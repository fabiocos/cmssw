/**
 * Unit test to compare ONNX vs TorchScript inference for dummy vertex model.
 * 
 * Build: scram b -j4
 * Run:   RecoVertex/PrimaryVertexProducer/test/testGNNBackends
 */

// Include PyTorch FIRST to avoid macro conflicts with Catch2
#include "PhysicsTools/PyTorch/interface/Model.h"

// Undef the CHECK macro from PyTorch before including Catch2
#ifdef CHECK
#undef CHECK
#endif

#include <catch2/catch_all.hpp>

#include <vector>
#include <cmath>
#include <iostream>

#include "FWCore/Utilities/interface/FileInPath.h"
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"

using cms::Ort::ONNXRuntime;

// Test input dimensions
constexpr int NUM_TRACKS = 50;
constexpr int NUM_FEATURES = 13;
constexpr int NUM_SLOTS = 180;  // Must match v17p1 model's num_slots

// Generate deterministic test input
std::vector<float> generateInput(int seed = 42) {
    std::vector<float> input(NUM_TRACKS * NUM_FEATURES);
    float val = static_cast<float>(seed);
    for (auto& x : input) {
        // Simple PRNG for reproducibility
        val = std::fmod(val * 1.1f + 0.3f, 10.0f);
        x = (val - 5.0f) / 5.0f;  // Range [-1, 1]
    }
    return input;
}

TEST_CASE("ONNX model loads and runs", "[gnn][onnx]") {
    const std::string modelPath = "RecoVertex/PrimaryVertexProducer/data/dummy_vertex_slot.onnx";
    
    SECTION("Model loads successfully") {
        REQUIRE_NOTHROW(ONNXRuntime(edm::FileInPath(modelPath).fullPath()));
    }
    
    SECTION("Inference produces correct output shapes") {
        ONNXRuntime onnx(edm::FileInPath(modelPath).fullPath());
        
        auto input = generateInput();
        std::vector<std::vector<float>> inputVec = {input};
        std::vector<std::vector<long int>> inputDims = {{1, NUM_TRACKS, NUM_FEATURES}};
        
        auto outputs = onnx.run({"x"}, inputVec, inputDims);
        
        // Check output count
        REQUIRE(outputs.size() == 5);
        
        // Check shapes: A[N,K], z_hat[K], t_hat[K], p[K], pi[K,4]
        CHECK(outputs[0].size() == static_cast<size_t>(NUM_TRACKS * NUM_SLOTS));  // A
        CHECK(outputs[1].size() == static_cast<size_t>(NUM_SLOTS));               // z_hat
        CHECK(outputs[2].size() == static_cast<size_t>(NUM_SLOTS));               // t_hat
        CHECK(outputs[3].size() == static_cast<size_t>(NUM_SLOTS));               // p
        CHECK(outputs[4].size() == static_cast<size_t>(NUM_SLOTS * 4));           // pi
    }
}

TEST_CASE("TorchScript model loads and runs", "[gnn][torchscript]") {
    const std::string modelPath = "RecoVertex/PrimaryVertexProducer/data/dummy_vertex_slot.pt";
    
    SECTION("Model loads successfully") {
        REQUIRE_NOTHROW(cms::torch::Model(edm::FileInPath(modelPath).fullPath()));
    }
    
    SECTION("Inference produces correct output shapes") {
        cms::torch::Model model(edm::FileInPath(modelPath).fullPath());
        
        // Create input tensor
        auto input = generateInput();
        auto options = ::torch::TensorOptions().dtype(::torch::kFloat32);
        auto inputTensor = ::torch::from_blob(input.data(), {1, NUM_TRACKS, NUM_FEATURES}, options);
        
        // Run inference
        std::vector<::torch::jit::IValue> inputs;
        inputs.push_back(inputTensor);
        auto output = model.forward(inputs);
        
        // Output is a tuple of 5 tensors
        REQUIRE(output.isTuple());
        auto tuple = output.toTuple();
        REQUIRE(tuple->elements().size() == 5);
        
        // Check shapes - use individual comparisons to avoid comma-in-macro issues
        auto A = tuple->elements()[0].toTensor();
        auto z_hat = tuple->elements()[1].toTensor();
        auto t_hat = tuple->elements()[2].toTensor();
        auto p = tuple->elements()[3].toTensor();
        auto pi = tuple->elements()[4].toTensor();
        
        CHECK(A.size(0) == NUM_TRACKS);
        CHECK(A.size(1) == NUM_SLOTS);
        CHECK(z_hat.size(0) == NUM_SLOTS);
        CHECK(t_hat.size(0) == NUM_SLOTS);
        CHECK(p.size(0) == NUM_SLOTS);
        CHECK(pi.size(0) == NUM_SLOTS);
        CHECK(pi.size(1) == 4);
    }
}

TEST_CASE("ONNX and TorchScript produce matching outputs", "[gnn][comparison]") {
    const std::string onnxPath = "RecoVertex/PrimaryVertexProducer/data/dummy_vertex_slot.onnx";
    const std::string tsPath = "RecoVertex/PrimaryVertexProducer/data/dummy_vertex_slot.pt";
    
    // Load models
    ONNXRuntime onnx(edm::FileInPath(onnxPath).fullPath());
    cms::torch::Model tsModel(edm::FileInPath(tsPath).fullPath());
    
    // Identical input
    auto input = generateInput(42);
    
    // ONNX inference
    std::vector<std::vector<float>> inputVec = {input};
    std::vector<std::vector<long int>> inputDims = {{1, NUM_TRACKS, NUM_FEATURES}};
    auto onnxOutputs = onnx.run({"x"}, inputVec, inputDims);
    
    // TorchScript inference
    auto options = ::torch::TensorOptions().dtype(::torch::kFloat32);
    auto inputTensor = ::torch::from_blob(input.data(), {1, NUM_TRACKS, NUM_FEATURES}, options);
    std::vector<::torch::jit::IValue> inputs;
    inputs.push_back(inputTensor);
    auto tsOutput = tsModel.forward(inputs);
    auto tuple = tsOutput.toTuple();
    
    // Compare each output
    const float tolerance = 1e-5f;
    
    SECTION("Assignment matrix A matches") {
        auto A_ts = tuple->elements()[0].toTensor().contiguous();
        auto A_ts_data = A_ts.data_ptr<float>();
        
        for (size_t i = 0; i < onnxOutputs[0].size(); ++i) {
            CHECK(std::abs(onnxOutputs[0][i] - A_ts_data[i]) < tolerance);
        }
    }
    
    SECTION("z_hat matches") {
        auto z_ts = tuple->elements()[1].toTensor().contiguous();
        auto z_ts_data = z_ts.data_ptr<float>();
        
        for (size_t i = 0; i < onnxOutputs[1].size(); ++i) {
            CHECK(std::abs(onnxOutputs[1][i] - z_ts_data[i]) < tolerance);
        }
    }
    
    SECTION("All outputs match within tolerance") {
        const char* names[] = {"A", "z_hat", "t_hat", "p", "pi"};
        bool all_match = true;
        
        for (int out_idx = 0; out_idx < 5; ++out_idx) {
            auto ts_tensor = tuple->elements()[out_idx].toTensor().contiguous();
            auto ts_data = ts_tensor.data_ptr<float>();
            
            float max_diff = 0;
            for (size_t i = 0; i < onnxOutputs[out_idx].size(); ++i) {
                max_diff = std::max(max_diff, std::abs(onnxOutputs[out_idx][i] - ts_data[i]));
            }
            
            std::cout << "  " << names[out_idx] << ": max_diff = " << max_diff << std::endl;
            if (max_diff >= tolerance) all_match = false;
        }
        
        CHECK(all_match);
    }
}
