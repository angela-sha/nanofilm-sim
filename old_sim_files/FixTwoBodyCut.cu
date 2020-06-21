#include "FixTwoBodyCut.h"

#include "BoundsGPU.h"
#include "GridGPU.h"
#include "list_macro.h"
#include "State.h"
#include "cutils_func.h"
#include "ReadConfig.h"
#include "EvaluatorWrapper.h"
#include "PairEvaluatorTwoBody.h"
#include "EvaluatorWrapper.h"
//#include "ChargeEvaluatorEwald.h"
using namespace std;
namespace py = boost::python;
const string TwoBodyCutType = "TwoBodyCut";



FixTwoBodyCut::FixTwoBodyCut(boost::shared_ptr<State> state_, string handle_)
    : FixPair(state_, handle_, "all", TwoBodyCutType, true, false, 1),
    epsHandle("eps"), sigHandle("sig"), rCutHandle("rCut"), sigma_RHandle("sigma_R"), r_AHandle("r_A"), alpha_AHandle("alpha_A"), GHandle("G"), r_GHandle(r_G), sigma_GHandle("sigma_G")
{

    initializeParameters(epsHandle, epsilons);
    initializeParameters(sigHandle, sigmas);
    initializeParameters(rCutHandle, rCuts);
    initializeParameters(sigma_RHandle, sigma_Rs);
    initializeParameters(r_AHandle, r_As);
    initializeParameters(alpha_AHandle, alpha_As);
    initializeParameters(GHandle, Gs);
    initializeParameters(r_GHandle, r_Gs);
    initializeParameters(sigma_GHandle, sigma_Gs);
    paramOrder = {rCutHandle, epsHandle, sigHandle, sigma_RHandle, r_AHandle, alpha_AHandle, GHandle, r_GHandle, sigma_GHandle};
    readFromRestart();
    canAcceptChargePairCalc = true;
    setEvalWrapper();
}

void FixTwoBodyCut::compute(int virialMode) {
    int nAtoms = state->atoms.size();
    int numTypes = state->atomParams.numTypes;
    GPUData &gpd = state->gpd;
    GridGPU &grid = state->gridGPU;
    int activeIdx = gpd.activeIdx();
    uint16_t *neighborCounts = grid.perAtomArray.d_data.data();
    float *neighborCoefs = state->specialNeighborCoefs;

    evalWrap->compute(nAtoms, gpd.xs(activeIdx), gpd.fs(activeIdx),
                      neighborCounts, grid.neighborlist.data(), grid.perBlockArray.d_data.data(),
                      state->devManager.prop.warpSize, paramsCoalesced.data(), numTypes, state->boundsGPU,
                      neighborCoefs[0], neighborCoefs[1], neighborCoefs[2], gpd.virials.d_data.data(), gpd.qs(activeIdx), chargeRCut, virialMode);

}

void FixTwoBodyCut::singlePointEng(float *perParticleEng) {
    int nAtoms = state->atoms.size();
    int numTypes = state->atomParams.numTypes;
    GPUData &gpd = state->gpd;
    GridGPU &grid = state->gridGPU;
    int activeIdx = gpd.activeIdx();
    uint16_t *neighborCounts = grid.perAtomArray.d_data.data();
    float *neighborCoefs = state->specialNeighborCoefs;
    evalWrap->energy(nAtoms, gpd.xs(activeIdx), perParticleEng, neighborCounts, grid.neighborlist.data(), grid.perBlockArray.d_data.data(), state->devManager.prop.warpSize, paramsCoalesced.data(), numTypes, state->boundsGPU, neighborCoefs[0], neighborCoefs[1], neighborCoefs[2], gpd.qs(activeIdx), chargeRCut);

}

void FixTwoBodyCut::setEvalWrapper() {
    if (evalWrapperMode == "offload") {
        EvaluatorTwoBody eval;
        evalWrap = pickEvaluator<EvaluatorTwoBody, 8, true>(eval, chargeCalcFix);
    } else if (evalWrapperMode == "self") {
        EvaluatorTwoBody eval;
        evalWrap = pickEvaluator<EvaluatorTwoBody, 8, true>(eval, nullptr);
    }
}

bool FixTwoBodyCut::prepareForRun() {
    //loop through all params and fill with appropriate lambda function, then send all to device
    auto fillEps = [] (float a, float b) {
        return sqrt(a*b);
    };

    auto fillSig = [] (float a, float b) {
            return (a+b) / 2.0;
    };
    auto fillRCut = [this] (float a, float b) {
        return (float) std::fmax(a, b);
    };
    auto none = [] (float a){};

    auto fillRCutDiag = [this] () {
        return (float) state->rCut;
    };

    auto processEps = [] (float a) {
        return 24*a;
    };
    auto processSig = [] (float a) {
        return pow(a, 6);
    };
    auto processRCut = [] (float a) {
        return a*a;
    };
    auto processSigma_R = [] (float a) {
        return a;
    };
    auto processR_A = [] (float a) {
        return a;
    };
    auto processAlpha_A = [] (float a) {
        return 1.0f / a;
    };
    auto processG = [] (float a) {
        return a;
    };
    auto processR_G = [] (float a) {
        return a;
    };
    auto processSigma_G = [] (float a) {
        return 2.0f / (a*a);
    };
    prepareParameters(epsHandle, fillEps, processEps, false);
    prepareParameters(sigHandle, fillSig, processSig, false);
    prepareParameters(rCutHandle, fillRCut, processRCut, true, fillRCutDiag);
    prepareParameters(sigma_RHandle, fillEps, processSigma_R, false);
    prepareParameters(r_AHandle, fillEps, processR_A, false);
    prepareParameters(alpha_AHandle, fillEps, processAlpha_A, false);
    prepareParameters(GHandle, fillEps, processG, false);
    prepareParameters(r_GHandle, fillEps, processR_G, false);
    prepareParameters(sigma_GHandle, fillEps, processSigma_G, false);

    sendAllToDevice();
    setEvalWrapper();
    return true;
}

string FixTwoBodyCut::restartChunk(string format) {
    stringstream ss;
    ss << restartChunkPairParams(format);
    return ss.str();
}


bool FixTwoBodyCut::postRun() {

    return true;
}

void FixTwoBodyCut::addSpecies(string handle) {
    initializeParameters(epsHandle, epsilons);
    initializeParameters(sigHandle, sigmas);
    initializeParameters(rCutHandle, rCuts);
    initializeParameters(sigma_RHandle, sigma_Rs);
    initializeParameters(r_AHandle, r_As);
    initializeParameters(alpha_AHandle, alpha_As);
    initializeParameters(GHandle, Gs);
    initializeParameters(r_GHandle, r_Gs);
    initializeParameters(sigma_GHandle, sigma_Gs);

}

vector<float> FixTwoBodyCut::getRCuts() {
    vector<float> res;
    vector<float> &src = *(paramMap[rCutHandle]);
    for (float x : src) {
        if (x == DEFAULT_FILL) {
            res.push_back(-1);
        } else {
            res.push_back(x);
        }
    }

    return res;
}

void export_FixTwoBodyCut() {
    py::class_<FixTwoBodyCut, boost::shared_ptr<FixTwoBodyCut>, py::bases<FixPair>, boost::noncopyable > (
        "FixTwoBodyCut",
        py::init<boost::shared_ptr<State>, string> (py::args("state", "handle"))
    )
      ;
}
