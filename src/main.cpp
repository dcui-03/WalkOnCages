#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"

#include <Eigen/Core>
#include <chrono>
#include <memory>
#include <iostream>
#include <string>
#include <cmath>
#include <cctype>
#include <algorithm>
#include <tuple>
#include <array>
#include <vector>
#include <stdexcept>
#include <glm/glm.hpp>
#include <glm/vec3.hpp>

//#include "args/args.hxx"
#include "imgui.h"

// My files
#include "cagedeformer/cagedeformer.hpp"
#include "cage/types/curvenetcage.hpp"
#include "cage/types/meshcage.hpp"
#include "query/types/meshquery.hpp"
#include "curvenet/curvenet.hpp"
#include "mesh/mesh.hpp"
#include "psCurvenet/pscurvenet.hpp"
#include "utils/utils.hpp"
#include "IO/io.hpp"

// Main file for visualization with Polyscope

// VARIABLES FOR POLYSCOPE OBJECTS

// SURFACE MESH (M)
// NOTE: this needs to be a generic surface mesh and NOT a triangle mesh
Eigen::MatrixXd psMesh_V; // Vertex list
std::vector<std::vector<int>> psMesh_F; // Face list: Note the inner list has arbitrary size for non-triangle faces
polyscope::SurfaceMesh* psMesh = nullptr;

// CAGE MESH topology (mesh-cage mode only). The cage is displayed as a wireframe + point
// cloud (via psEditableCN/psControlsPC below) rather than a surface mesh, so it doesn't
// visually cover the query mesh.
std::vector<std::vector<int>> psCageMesh_F;
std::vector<std::array<int, 2>> psCageWireframe_E;

// CURVENET (CN) / CAGE WIREFRAME (mesh-cage mode)
Eigen::MatrixXd psCN_P; // Point list
std::vector<std::array<int, 2>> psCN_E; // Edge List
polyscope::CurveNetwork* psEditableCN = nullptr;

// CONTROLS (PC) / CAGE POINTS (mesh-cage mode). This is the point cloud the user drags in
// either mode: curvenet controls, or cage mesh vertices.
Eigen::MatrixXd psControls_P; // Controls
polyscope::PointCloud* psControlsPC = nullptr;    // point cloud for controls
polyscope::PointCloudColorQuantity* psCageColorQ = nullptr; // cage's own colors, for comparison
// TANGENTS (PC) - curvenet mode only, no equivalent for mesh cages
Eigen::MatrixXd psTangents_P; // Tangents
polyscope::PointCloud* psTangentsPC = nullptr;    // point cloud for tangents

// TANGENTS (CN)
Eigen::MatrixXd psTangentsVec; // Aggregate list of controls and tangents
std::vector<std::array<int, 2>> psTangents_E; // Edge List between controls and tangents
polyscope::CurveNetwork* psTangentsCN = nullptr;  // Connects controls to their tangents


// VARIABLES FOR PARSING AND WRITING FILES
std::string InputPath;
std::string OutputPath;
std::string CurvenetPath;
std::string CageMeshPath;
bool loadedCurvenet = false;
// True if a cage mesh OBJ was given on the command line; false means curvenet-cage mode
bool meshCageMode = false;

// UI HELPERS
bool CD_init = false;

bool disable_psCN = false;

// Real-time dragging updates either positions or debug colors/gradients, never both at once
bool colorMode = false;
bool debugColorsVisible = true;

bool createCtrlMode = false;  // Allows users to place control points
bool createSplineMode = false;  // Allow users to initialize new splines

bool editCtrlMode = false;   // Allows users to modify controls / cage points
bool editTanMode = false;   // Allows users to modify tangents (curvenet mode only)

bool delCtrlMode = false;   // Allows users to remove control points
bool delSplineMode = false;  // Allows users to remove splines

// Spline creation/removal helpers
int selectedIdx = -1;    // Index of selected vertex on mesh
std::pair<int, int> selectedPair = {-1, -1};

// Editing helpers
bool tanConstraint = true;  // Constrain tangent movement to tangent plane only

// Gizmo helpers
bool activeGizmo = false; // This tells us if there is an active gizmo
Eigen::Vector3d gizmoPos;
static polyscope::TransformationGizmo* vertexGizmo = nullptr;

// Pre-computation
int samplingParam = 5;    // Curvenet discretization density (alpha)
int num_samples = 50;     // WoS samples per query point
float offsetParam = 0.5f; // Surface offset applied to newly-created controls

// Discrete curvenet for modeling
std::unique_ptr<psCurvenet::pscurvenet> psCN = nullptr; // Curvenet that polyscope will use for updates

// Cage Deformer
std::unique_ptr<Mesh::mesh> CD_Mesh = nullptr;         // Query mesh
std::unique_ptr<Mesh::mesh> CageMesh = nullptr;        // Cage mesh (mesh-cage mode only)
std::unique_ptr<Curvenet::curvenet> CD_CN = nullptr;   // Curvenet-cage mode only
std::unique_ptr<Cage::cage> CD_Cage = nullptr;         // Either a curvenetcage or a meshcage
std::unique_ptr<Query::meshquery> CD_Query = nullptr;
std::unique_ptr<CageDeformer::cagedeformer> CD = nullptr;

// Debug color/gradient quantity handles; nulled in resetMesh() when the structure is dropped
polyscope::SurfaceVertexColorQuantity* psColorQ = nullptr;
polyscope::SurfaceVertexVectorQuantity* psGradRQ = nullptr;
polyscope::SurfaceVertexVectorQuantity* psGradGQ = nullptr;
polyscope::SurfaceVertexVectorQuantity* psGradBQ = nullptr;


// ----------------- FUNCTIONS BEGIN HERE -------------------------

// Saves current curvenet to some file format
int saveCurvenet() {
    if (!psCN) {
        std::cout << "No curvenet object to save." << std::endl;
        return -1;
    }

    if (OutputPath.empty()) {
        std::cout << "No output .curvenet path specified." << std::endl;
        return -1;
    }

    int success = psCN->saveCurvenet(OutputPath);

    if (success == 1) {
        std::cout << "Curvenet saved to: " << OutputPath << std::endl;
    } else {
        std::cout << "Failed to save curvenet to: " << OutputPath << std::endl;
    }

    return success;
}

// Prints a rejection message and returns true if the caller should stop, since this
// action needs an actual curvenet and we're in mesh-cage mode
bool rejectMeshCage(const std::string& reason) {
    if (!meshCageMode) {
        return false;
    }
    std::cout << "Using mesh cage mode. " << reason << std::endl;
    return true;
}

// Creates gizmo at vertex
void addGizmoAtLocation(Eigen::Vector3d& startpos) {
    activeGizmo = true;
    if (!vertexGizmo) {
        vertexGizmo = polyscope::addTransformationGizmo("vertex_editor");
        vertexGizmo->setAllowTranslation(true);
        vertexGizmo->setAllowRotation(false);
        vertexGizmo->setAllowScaling(false);
        vertexGizmo->setInteractInLocalSpace(false);
        //vertexGizmo->setGizmoSize(0.5f);
    }

    // Place gizmo at vertex position
    vertexGizmo->setPosition(Utils::eigenToGLM(startpos));
    return;
}
// Creates a new gizmo rotated to match an axis
void addGizmoAtLocation(Eigen::Vector3d& startpos, Eigen::Vector3d& ax1, Eigen::Vector3d& ax2, Eigen::Vector3d& ax3) {
    activeGizmo = true;
    if (!vertexGizmo) {
        glm::mat4 T(1.0f);
        T[0] = glm::vec4(Utils::eigenToGLM(ax1), 0.0f);
        T[1] = glm::vec4(Utils::eigenToGLM(ax2), 0.0f);
        T[2] = glm::vec4(Utils::eigenToGLM(ax3), 0.0f);
        T[3] = glm::vec4(Utils::eigenToGLM(startpos), 1.0f);
        vertexGizmo = polyscope::addTransformationGizmo("editor");
        vertexGizmo->setTransform(T);
        vertexGizmo->setAllowTranslation(true);
        vertexGizmo->setAllowRotation(true);
        vertexGizmo->setAllowScaling(false);
        vertexGizmo->setInteractInLocalSpace(true); // TODO: Check this
        //vertexGizmo->setGizmoSize(0.5f);
    }

    // place gizmo at vertex position
    vertexGizmo->setPosition(Utils::eigenToGLM(startpos));
    gizmoPos = startpos;
    return;
}

// Removes gizmo
void removeGizmo() {
    if (!activeGizmo) return;

    if (vertexGizmo) {
        vertexGizmo->remove();
        vertexGizmo = nullptr;
    }
    activeGizmo = false;
    return;
}

// Remove all polyscope object
void removeAllCurvenetPS() {
    if (psEditableCN) {
        psEditableCN->remove();
        psEditableCN = nullptr;
    }

    if (psTangentsCN) {
        psTangentsCN->remove();
        psTangentsCN = nullptr;
    }

    if (psControlsPC) {
        psControlsPC->remove();
        psControlsPC = nullptr;
        psCageColorQ = nullptr;
    }

    if (psTangentsPC) {
        psTangentsPC->remove();
        psTangentsPC = nullptr;
    }
    return;
}

// Update the mesh vertex positions
void updateMesh(const std::vector<Eigen::Vector3d>& new_pos) {
    // Convert to matrix form
    if (new_pos.size() != psMesh_V.rows()) {
        std::cout << "Invalid mesh update size: New pos has size " << new_pos.size() << " but mesh has size " << psMesh_V.rows() << std::endl;
        return;
    }
    for (int v = 0; v < psMesh_V.rows(); v++) {
        psMesh_V.row(v) = new_pos[v].transpose();
    }
    // Update mesh
    psMesh->updateVertexPositions(psMesh_V);
    return;
}

void setDebugColorsVisible(bool visible) {
    debugColorsVisible = visible;
    if (psColorQ) psColorQ->setEnabled(debugColorsVisible);
    if (psGradRQ) psGradRQ->setEnabled(debugColorsVisible);
    if (psGradGQ) psGradGQ->setEnabled(debugColorsVisible);
    if (psGradBQ) psGradBQ->setEnabled(debugColorsVisible);
    if (psCageColorQ) psCageColorQ->setEnabled(debugColorsVisible);
    return;
}

// Deformation mode: recompute mesh positions
void updateCageDeformerPositions(bool recompute = true) {
    if (!CD_init || !CD) {
        return;
    }
    if (recompute) {
        if (!meshCageMode) {
            std::vector<Eigen::Vector3d> controlsV, tangentsV;
            std::vector<std::array<int, 4>> splines;
            psCN->cnAsStdVector(controlsV, tangentsV, splines);
            CD_CN->updateCurveNet(controlsV, tangentsV);
        }

        Eigen::MatrixXd new_pos;
        if (CD->applyDeformation(new_pos) == 1) {
            std::vector<Eigen::Vector3d> new_pos_v;
            Utils::EigM3toStdV(new_pos, new_pos_v);
            updateMesh(new_pos_v);
        }
    }
    return;
}

// Color mode: recompute debug colors/gradients only
void updateCageDeformerColors(bool recompute = true) {
    if (!CD_init || !CD) {
        return;
    }
    if (recompute) {
        if (!meshCageMode) {
            std::vector<Eigen::Vector3d> controlsV, tangentsV;
            std::vector<std::array<int, 4>> splines;
            psCN->cnAsStdVector(controlsV, tangentsV, splines);
            CD_CN->updateCurveNet(controlsV, tangentsV);
        }

        std::vector<glm::vec3> colors;
        if (CD->colorsPolyscopeFormat(colors) == 1) {
            psColorQ = psMesh->addVertexColorQuantity("Stochastic Colors", colors);
        } else {
            std::cout << "Failed to compute debug colors." << std::endl;
        }
        std::vector<glm::vec3> gradR, gradG, gradB;
        if (CD->colorGradientsPolyscopeFormat(gradR, gradG, gradB) == 1) {
            psGradRQ = psMesh->addVertexVectorQuantity("Color Gradient R", gradR);
            psGradRQ->setVectorColor(glm::vec3{1.0f, 0.0f, 0.0f});
            psGradGQ = psMesh->addVertexVectorQuantity("Color Gradient G", gradG);
            psGradGQ->setVectorColor(glm::vec3{0.0f, 1.0f, 0.0f});
            psGradBQ = psMesh->addVertexVectorQuantity("Color Gradient B", gradB);
            psGradBQ->setVectorColor(glm::vec3{0.0f, 0.0f, 1.0f});
        } else {
            std::cout << "Failed to compute debug color gradients." << std::endl;
        }
        // Cage's own colors, for comparison. Curvenet cages give colors for controls+tangents;
        // the point cloud only shows controls, so take the leading slice
        std::vector<glm::vec3> cageColors;
        int numPts = static_cast<int>(psControls_P.rows());
        if (psControlsPC && CD->cageColorsPolyscopeFormat(cageColors) == 1 && static_cast<int>(cageColors.size()) >= numPts) {
            std::vector<glm::vec3> ptColors(cageColors.begin(), cageColors.begin() + numPts);
            psCageColorQ = psControlsPC->addColorQuantity("Cage Colors", ptColors);
        } else {
            std::cout << "Failed to compute cage debug colors." << std::endl;
        }
        // Re-registering may reset enabled state; re-apply it
        setDebugColorsVisible(debugColorsVisible);
    }
    return;
}

// Update the curvenet object
void updateCurvenet(bool conn = false) {
    // Reset curvenet
    psCN->cnAsCurveNetwork(psCN_P, psCN_E);
    psCN->tansAsCurveNetwork(psTangentsVec, psTangents_E);
    psCN->cPosAsMatrix(psControls_P);
    psCN->tPosAsMatrix(psTangents_P);

    if (conn) {
        removeAllCurvenetPS();
        if (psCN_E.size() > 0) {
            psEditableCN = polyscope::registerCurveNetwork("Curvenet", psCN_P, psCN_E);
            psEditableCN->setColor({0.0f, 0.0f, 1.0f});
            psEditableCN->setMaterial("flat");
            psEditableCN->setTransparency(0.65);
            psEditableCN->setRadius(0.003);
            psEditableCN->setEnabled(true);
        }

        if (psTangents_E.size() > 0) {
            psTangentsCN = polyscope::registerCurveNetwork("Handles", psTangentsVec, psTangents_E);
            psTangentsCN->setColor({0.5f, 0.55f, 0.15f});
            psTangentsCN->setMaterial("flat");
            psTangentsCN->setTransparency(0.8);
            psTangentsCN->setRadius(0.004);
            psTangentsCN->setEnabled(true);
        }

        if (psControls_P.rows() > 0) {
            psControlsPC = polyscope::registerPointCloud("Controls", psControls_P);
            psControlsPC->setPointColor({0.9f, 0.2f, 0.1f});
            psControlsPC->setMaterial("flat");
            psControlsPC->setPointRadius(0.02);
            psControlsPC->setEnabled(true);
        }

        if (psTangents_P.rows() > 0) {
            psTangentsPC = polyscope::registerPointCloud("Tangents", psTangents_P);
            psTangentsPC->setPointColor({0.1f, 0.9f, 0.2f});
            psTangentsPC->setMaterial("flat");
            psTangentsPC->setPointRadius(0.012);
            psTangentsPC->setEnabled(true);
        }
    } else {
        if (psEditableCN && psCN_E.size() > 0) {
            psEditableCN->updateNodePositions(psCN_P);
        }
        if (psTangentsCN && psTangents_E.size() > 0) {
            psTangentsCN->updateNodePositions(psTangentsVec);
        }
        if (psControlsPC && psControls_P.rows() > 0) {
            psControlsPC->updatePointPositions(psControls_P);
        }
        if (psTangentsPC && psTangents_P.rows() > 0) {
            psTangentsPC->updatePointPositions(psTangents_P);
        }
    }
    return;
}

// Mesh-cage mode's analog of updateCurvenet: shows the cage as a wireframe (psEditableCN)
// + point cloud (psControlsPC), sharing the same polyscope objects curvenet mode uses.
// Topology is fixed (loaded once from OBJ), so this only ever moves positions.
void updateCageMeshViz(bool conn = false) {
    if (conn) {
        removeAllCurvenetPS();
        psEditableCN = polyscope::registerCurveNetwork("Cage Wireframe", psControls_P, psCageWireframe_E);
        psEditableCN->setColor({0.0f, 0.0f, 1.0f});
        psEditableCN->setMaterial("flat");
        psEditableCN->setTransparency(0.65);
        psEditableCN->setRadius(0.003);
        psEditableCN->setEnabled(true);

        psControlsPC = polyscope::registerPointCloud("Cage Points", psControls_P);
        psControlsPC->setPointColor({0.9f, 0.2f, 0.1f});
        psControlsPC->setMaterial("flat");
        psControlsPC->setPointRadius(0.02);
        psControlsPC->setEnabled(true);
    } else {
        if (psEditableCN) {
            psEditableCN->updateNodePositions(psControls_P);
        }
        if (psControlsPC) {
            psControlsPC->updatePointPositions(psControls_P);
        }
    }
    return;
}

void resetMesh() {
    if (!IO::readOBJ(InputPath, psMesh_V, psMesh_F)) {
        return;
    }
    psMesh = polyscope::registerSurfaceMesh("Surface Mesh", psMesh_V, psMesh_F);
    psMesh->setSurfaceColor({0.6f, 0.6f, 0.6f});
    // Old quantity handles would dangle otherwise
    psColorQ = nullptr;
    psGradRQ = nullptr;
    psGradGQ = nullptr;
    psGradBQ = nullptr;
    return;
}

void clearCD() {
    CD_init = false;
    CD = nullptr;
    CD_Cage = nullptr;
    CD_Query = nullptr;
    CD_CN = nullptr;
    return;
}

// Shared tail of "Compute CageDeformer": assumes CD_Cage is already set
void finalizeCageDeformer() {
    CD_Query = std::make_unique<Query::meshquery>(CD_Mesh.get());
    CD = std::make_unique<CageDeformer::cagedeformer>();
    CD->applyCage(CD_Cage.get());
    CD->applyQuery(CD_Query.get());
    int success = CD->computeCoordinates(0, num_samples);
    if (success != 1) {
        std::cout << "Failed to compute stochastic barycentric coordinates." << std::endl;
        clearCD();
    } else {
        CD_init = true;
        updateCageDeformerColors(true);
    }
    return;
}

void disableCurvenet() {
    if (psControls_P.rows() < 1) {
        disable_psCN = true;
        std::cout << "No curvenet object to enable/disable." << std::endl;
        return;
    }

    if (disable_psCN) {
        if (psTangentsCN) {
            psTangentsCN->setEnabled(false);
        }
        if (psTangentsPC) {
            psTangentsPC->setEnabled(false);
        }
        if (psControlsPC) {
            psControlsPC->setEnabled(false);
        }
        if (psEditableCN) {
            psEditableCN->setEnabled(false);
        }
        disable_psCN = false;
    } else {
        if (psTangentsCN) {
            psTangentsCN->setEnabled(true);
        }
        if (psTangentsPC) {
            psTangentsPC->setEnabled(true);
        }
        if (psControlsPC) {
            psControlsPC->setEnabled(true);
        }
        if (psEditableCN) {
            psEditableCN->setEnabled(true);
        }
        disable_psCN = true;
    }
}
// clear all modes and their variables except for the specified mode
int clearModes() {
    createCtrlMode = false;
    createSplineMode = false;
    editCtrlMode = false;
    editTanMode = false;
    delCtrlMode = false;
    delSplineMode = false;

    removeGizmo();

    selectedIdx = -1;
    gizmoPos = Eigen::Vector3d::Zero();

    selectedPair = {-1, -1};
    return 1;
}

void ImGuiSection(const char* label) {
    ImGui::Spacing();
    ImGui::Separator();
    ImGui::TextUnformatted(label);
    ImGui::Spacing();
}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {
    ImGuiIO& io = ImGui::GetIO();
    bool mouseDown = ImGui::IsMouseDown(0);
    bool mouseClicked = ImGui::IsMouseClicked(0);
    bool mouseReleased = ImGui::IsMouseReleased(0);
    glm::vec2 screen{io.MousePos.x, io.MousePos.y};

    polyscope::PickResult pick = polyscope::pickAtScreenCoords(screen);

    ImGuiSection("Cage Deformer Initialization");
    // Pre-compute the cage, query, and stochastic barycentric coordinates
    if (ImGui::Button("Compute CageDeformer")) {
        clearModes();
        if (meshCageMode) {
            clearCD();
            // Rebuild fresh from the current (possibly drag-edited) point cloud, so the
            // BVH and all derived mesh data reflect the latest positions
            std::vector<Eigen::Vector3d> cageV;
            Utils::EigM3toStdV(psControls_P, cageV);
            CageMesh = std::make_unique<Mesh::mesh>(cageV, psCageMesh_F);
            CD_Cage = std::make_unique<Cage::meshcage>(CageMesh.get());
            finalizeCageDeformer();
        } else if (psControls_P.rows() <= 1) {
            std::cout << "No splines specified. Add one or more splines." << std::endl;
        } else {
            // Compress the curvenet
            psCN->cleanupControls();
            saveCurvenet();
            clearCD();
            updateCurvenet();
            std::vector<Eigen::Vector3d> controlsV, tangentsV;
            std::vector<std::array<int, 4>> splines;
            psCN->cnAsStdVector(controlsV, tangentsV, splines);
            CD_CN = std::make_unique<Curvenet::curvenet>(controlsV, tangentsV, splines, CD_Mesh.get(), samplingParam);
            CD_Cage = std::make_unique<Cage::curvenetcage>(CD_CN.get());
            finalizeCageDeformer();
        }
    }

    // User parameters
    ImGui::SliderInt("Sampling Param", &samplingParam, 2, 8);
    ImGui::SliderInt("Num Samples", &num_samples, 5, 50);
    num_samples = std::clamp(((num_samples + 2) / 5) * 5, 5, 50);

    if (ImGui::Button(colorMode ? "Switch to Deformation Mode" : "Switch to Color Mode")) {
        colorMode = !colorMode;
    }
    ImGui::SameLine();
    ImGui::TextUnformatted(colorMode ? "(dragging updates: colors)" : "(dragging updates: positions)");
    if (ImGui::Button(debugColorsVisible ? "Hide Debug Colors" : "Show Debug Colors")) {
        setDebugColorsVisible(!debugColorsVisible);
    }

    ImGuiSection("Saving and Viewing");
    if (!meshCageMode) {
        // Save the current curvenet state
        if (ImGui::Button("Save Curvenet")) {
            clearModes();
            saveCurvenet();
        }
        ImGui::SameLine();
    }
    if (ImGui::Button("Save Mesh")) {   // TODO
        clearModes();
    }
    if (!meshCageMode) {
        if (ImGui::Button(disable_psCN ? "Disable Curvenet" : "Enable Curvenet")) {
            disableCurvenet();
        }
    }

    ImGuiSection("Create and Edit Splines");
    ImGui::SliderFloat("Offset", &offsetParam, 0.1f, 1.0f);
    // CONTROL/SPLINE CREATION
    // Create controls
    if (ImGui::Button(createCtrlMode ? "Stop Creating Controls" : "Create Controls")) {
        if (!rejectMeshCage("Cannot create curvenet.")) {
            bool tempMode = createCtrlMode;
            clearModes();
            createCtrlMode = !tempMode;
        }
    }
    ImGui::SameLine();
    // Move a curvenet vertex
    if (ImGui::Button(createSplineMode ? "Stop Creating Splines" : "Create Splines")) {
        if (!rejectMeshCage("Cannot create curvenet.")) {
            bool tempMode = createSplineMode;
            clearModes();
            createSplineMode = !tempMode;
        }
    }


    // CONTROL/TANGENT EDITING (Move Control is shared between curvenet controls and cage points)
    if (ImGui::Button(editCtrlMode ? "Stop Moving Control" : "Move Control")) {
        bool tempMode = editCtrlMode;
        clearModes();
        if (psControls_P.rows() == 0) {
            std::cout << "Cannot edit. No existing controls." << std::endl;
            editCtrlMode = false;
        } else if (tempMode == false) {
            editCtrlMode = true;
        }
    }
    ImGui::SameLine();
    if (ImGui::Button(editTanMode ? "Stop Moving Splines" : "Move Splines")) {
        if (!rejectMeshCage("No tangents to modify.")) {
            bool tempMode = editTanMode;
            clearModes();
            if (psTangents_P.rows() == 0) {
                std::cout << "Cannot edit. No existing splines." << std::endl;
                editTanMode = !tempMode;
            } else if (tempMode == false) {
                editTanMode = true;
            }
        }
    }
    ImGui::SameLine();
    ImGui::Checkbox("Proj. Tans", &tanConstraint); // TODO: do not allow degenerate vectors --> Constrain tan vertex AND gizmo

    // CONTROL/SPLINE DELETION
    if (ImGui::Button(delCtrlMode ? "Stop Removing Controls" : "Remove Controls")) {
        if (!rejectMeshCage("Cannot create curvenet.")) {
            bool tempMode = delCtrlMode;
            clearModes();
            if (psControls_P.rows() == 0) {
                std::cout << "Cannot delete. No existing controls." << std::endl;
                delCtrlMode = false;
            } else if (tempMode == false) {
                delCtrlMode = true;
            }
        }
    }
    ImGui::SameLine();
    if (ImGui::Button(delSplineMode ? "Stop Removing Splines" : "Remove Splines")) {
        if (!rejectMeshCage("Cannot create curvenet.")) {
            bool tempMode = delSplineMode;
            clearModes();
            if (psTangents_P.rows() == 0) {
                std::cout << "Cannot delete. No existing splines." << std::endl;
                delSplineMode = !tempMode;
            } else if (tempMode == false) {
                delSplineMode = true;
            }
        }
    }


    ImGuiSection("Resets");
    // RESETs
    if (ImGui::Button("Clear Gizmo")) {
        std::cout << "Removing current gizmo." << std::endl;
        removeGizmo();
        selectedIdx = -1;
    }
    if (!meshCageMode) {
        ImGui::SameLine();
        if (ImGui::Button("Clear Curvenet")) {
            std::cout << "Clearing entire curvenet." << std::endl;
            psCN->resetCurvenet();
            std::cout << "Clearing cage deformer object." << std::endl;
            clearCD();
            updateCurvenet(true);
            clearModes();
            // Reset mesh
            std::cout << "Resetting mesh." << std::endl;
            resetMesh();
            psMesh->setEnabled(true);
        }
    }
    // May need to store a copy of the rest curvenet
    if (ImGui::Button("Clear CageDeformer")) {
        std::cout << "Clearing cage deformer." << std::endl;
        clearCD();
        resetMesh();
        clearModes();
    }

    // CREATE MODE CLICKS
    // Clicked on mesh during create mode
    if (createCtrlMode && mouseClicked && pick.isHit) {
        // If we click a current control, intializes the control and tangents
        if (pick.structure == psMesh) {
            polyscope::SurfaceMeshPickResult meshPick = psMesh->interpretPickResult(pick);
            Eigen::Vector3d pos = Utils::glmToEigen(pick.position);

            // TODO: Switch to a different version for non-triangle meshes
            Eigen::Vector3d proj;
            Mesh::vertProjData vProjData = CD_Mesh->computeVProjection(pos, proj);
            if (vProjData.elIdx != -1 && vProjData.elType != -1) {
                Eigen::Vector3d normal = CD_Mesh->getNormal(vProjData);
                psCN->addControl(pos, normal, static_cast<double>(offsetParam));
                clearCD();
                updateCurvenet(true);
                std::cout << "New Vert created at (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << std::endl;
            } else {
                std::cout << "No valid point picked." << std::endl;
            }
        }
    }
    // Create Spline by picking two controls
    if (createSplineMode && mouseClicked && pick.isHit) {
        if (pick.structure == psControlsPC) {
            polyscope::PointCloudPickResult pcPick = psControlsPC->interpretPickResult(pick);

            selectedIdx = static_cast<int>(pcPick.index);

            if (selectedPair.first == -1) {
                selectedPair.first = selectedIdx;
                Eigen::Vector3d pos = psControls_P.row(selectedIdx).transpose();
                std::cout << "First Spline Vert: (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << std::endl;
            } else {
                selectedPair.second = selectedIdx;
                Eigen::Vector3d pos = psControls_P.row(selectedIdx).transpose();
                std::cout << "Second Spline Vert: (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << std::endl;
                // Compute initial tangent directions
                psCN->addSpline(selectedPair.first, selectedPair.second);
                // Reset pair
                selectedPair = {-1, -1};
                std::cout << "New Spline Created.\n" << std::endl;
                clearCD();
                updateCurvenet(true);
            }
        }
    }

    // REMOVE MODE CLICKS
    // Clicked on control to remove
    if (delCtrlMode && mouseClicked && pick.isHit && pick.structure == psControlsPC) {
        polyscope::PointCloudPickResult pcPick = psControlsPC->interpretPickResult(pick);

        psCN->removeControl(static_cast<int>(pcPick.index));
        std::cout << "Control removed." << std::endl;
        clearCD();
        updateCurvenet(true);
    }
    // Clicked on a tangent whose spline we should remove
    if (delSplineMode && mouseClicked && pick.isHit && pick.structure == psTangentsPC) {
        polyscope::PointCloudPickResult pcPick = psTangentsPC->interpretPickResult(pick);

        psCN->removeSplineByTangent(static_cast<int>(pcPick.index));
        std::cout << "Spline removed." << std::endl;
        clearCD();
        updateCurvenet(true);
    }

    // EDIT MODE CLICKS
    // Select control / cage point to edit
    if (editCtrlMode && mouseClicked) {
        if (pick.isHit && pick.structure == psControlsPC) {
            polyscope::PointCloudPickResult pcPick = psControlsPC->interpretPickResult(pick);

            selectedIdx = static_cast<int>(pcPick.index);
            // Build a basis
            Eigen::Vector3d selectedPos = psControls_P.row(selectedIdx).transpose();
            Eigen::Vector3d selectedN = meshCageMode ? CageMesh->getVNormal(selectedIdx) : psCN->getNormal(selectedIdx);
            Eigen::Vector3d t0, t1;
            Utils::buildPlaneBasis(selectedN, t0, t1);
            std::cout << "Editing Vert at (" << selectedPos[0] << ", " << selectedPos[1] << ", " << selectedPos[2] << ")" << std::endl;
            // add Gizmo at position
            addGizmoAtLocation(selectedPos, selectedN, t0, t1);
        } else if (editCtrlMode && mouseClicked && !activeGizmo && (!pick.isHit || (pick.isHit && pick.structure != psControlsPC))) {  // or clear
            clearModes();
            editCtrlMode = true;
        }
    }
    // Select tangent to edit (curvenet mode only; button rejects in mesh-cage mode)
    if (editTanMode && mouseClicked) {
        if (pick.isHit && pick.structure == psTangentsPC) {
            // Index into tangent list. Get associated spline by integer dividing by 2.
            // Then use the spline index to edit the tangent's position directly
            polyscope::PointCloudPickResult pcPick = psTangentsPC->interpretPickResult(pick);

            selectedIdx = static_cast<int>(pcPick.index);
            Eigen::Vector3d selectedPos = psTangents_P.row(selectedIdx).transpose();
            std::cout << "Editing Spline Handle at (" << selectedPos[0] << ", " << selectedPos[1] << ", " << selectedPos[2] << ")" << std::endl;
            // add Gizmo at position
            addGizmoAtLocation(selectedPos);
        } else if (editTanMode && mouseClicked && !activeGizmo && (!pick.isHit || (pick.isHit && pick.structure != psTangentsPC))) {  // or clear
            clearModes();
            editTanMode = true;
        }
    }

    // Update control / cage point position per-frame
    if (editCtrlMode && activeGizmo && selectedIdx >= 0) {
        // Get the gizmo's location at this frame
        Eigen::Vector3d gizmoPosF = Utils::glmToEigen(vertexGizmo->getPosition());

        if (meshCageMode) {
            psControls_P.row(selectedIdx) = gizmoPosF.transpose();
            CageMesh->setVertPos(selectedIdx, gizmoPosF);
            vertexGizmo->setPosition(Utils::eigenToGLM(gizmoPosF));
            updateCageMeshViz();
        } else {
            glm::mat4 T = vertexGizmo->getTransform();
            Eigen::Vector3d gizmoNormal = Utils::glmToEigen(glm::normalize(glm::vec3(T[0])));
            psCN->updateControlPos(selectedIdx, gizmoPosF, false);
            psCN->updateControlNormal(selectedIdx, gizmoNormal, true, false);
            vertexGizmo->setPosition(Utils::eigenToGLM(gizmoPosF));
            updateCurvenet();
        }

        if (colorMode) {
            updateCageDeformerColors(true);
        } else {
            updateCageDeformerPositions(true);
        }
    }
    // Update tangent position per-frame (curvenet mode only)
    if (editTanMode && activeGizmo && selectedIdx >= 0) {
        // Get the gizmo's location at this frame
        Eigen::Vector3d gizmoPosF = Utils::glmToEigen(vertexGizmo->getPosition());
        // If we are too close to either endpoint, do not update
        bool updated = psCN->updateTangentPos(selectedIdx, gizmoPosF, tanConstraint);
        updateCurvenet();
        if (colorMode) {
            updateCageDeformerColors(true);
        } else {
            updateCageDeformerPositions(true);
        }
        if (updated) {
            Eigen::Vector3d tangentPos = psTangents_P.row(selectedIdx).transpose();
            vertexGizmo->setPosition(Utils::eigenToGLM(tangentPos));
        }
    }

    return;
}

// Case-insensitive check for a path's extension
bool hasExtension(const std::string& path, const std::string& ext) {
    if (path.size() < ext.size()) {
        return false;
    }
    std::string suffix = path.substr(path.size() - ext.size());
    std::transform(suffix.begin(), suffix.end(), suffix.begin(), [](unsigned char c) { return std::tolower(c); });
    return suffix == ext;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        std::cout << "Too few arguments.\n"
                  << "Usage: ./stochastic_bc <query mesh OBJ> <cage mesh OBJ>\n"
                  << "   or: ./stochastic_bc <query mesh OBJ> <curvenet output path> [--load <curvenet path>]"
                  << std::endl;
        return 1;
    }
    InputPath = argv[1];
    std::string arg2 = argv[2];
    meshCageMode = hasExtension(arg2, ".obj");
    if (meshCageMode) {
        CageMeshPath = arg2;
    } else {
        OutputPath = arg2;
    }
    for (int i = 3; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "--load") {
            if (i + 1 >= argc) {
                std::cout << "Missing path after --load." << std::endl;
                return 1;
            }
            if (meshCageMode) {
                std::cout << "--load is ignored when a cage mesh is provided." << std::endl;
                i++;
            } else {
                CurvenetPath = argv[++i];
                loadedCurvenet = true;
            }
        } else {
            std::cout << "Unknown argument: " << arg << std::endl;
            return 1;
        }
    }

    // Initialize polyscope
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None; // Disable ground plane

    // Set the callback function
    polyscope::state::userCallback = myCallback;

    polyscope::init();

    // Load the query mesh
    std::cout << "\nLoading query mesh file" << std::endl;
    if (!IO::readOBJ(InputPath, psMesh_V, psMesh_F)) {
        std::cout << "Could not read query mesh" << std::endl;
        return -1;
    }
    psMesh = polyscope::registerSurfaceMesh("Surface Mesh", psMesh_V, psMesh_F);
    psMesh->setSurfaceColor({0.6f, 0.6f, 0.6f});

    std::vector<Eigen::Vector3d> meshV;
    Utils::EigM3toStdV(psMesh_V, meshV);
    CD_Mesh = std::make_unique<Mesh::mesh>(meshV, psMesh_F);

    psCN = std::make_unique<psCurvenet::pscurvenet>();

    if (meshCageMode) {
        std::cout << "Loading cage mesh file" << std::endl;
        if (!IO::readOBJ(CageMeshPath, psControls_P, psCageMesh_F)) {
            std::cout << "Could not read cage mesh" << std::endl;
            return -1;
        }
        IO::facesToWireframe(psCageMesh_F, psCageWireframe_E);
        updateCageMeshViz(true);

        std::vector<Eigen::Vector3d> cageMeshV;
        Utils::EigM3toStdV(psControls_P, cageMeshV);
        CageMesh = std::make_unique<Mesh::mesh>(cageMeshV, psCageMesh_F);
    } else if (loadedCurvenet) {
        int loadSuccess = psCN->loadCurvenet(CurvenetPath);
        if (loadSuccess == 1) {
            std::cout << "Loaded curvenet from: " << CurvenetPath << std::endl;
            updateCurvenet(true);
        } else {
            std::cout << "Failed to load curvenet from: " << CurvenetPath
                    << ". Starting from an empty curvenet." << std::endl;
            psCN->resetCurvenet();
        }
    }

    // Give control to the polyscope gui
    polyscope::show();

    return EXIT_SUCCESS;
}
