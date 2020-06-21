#include "FixRLCut.h"

#include "BoundsGPU.h"
#include "GridGPU.h"
#include "list_macro.h"
#include "State.h"
#include "cutils_func.h"
#include "ReadConfig.h"
#include "EvaluatorWrapper.h"
#include "PairEvaluatorRL.h"
#include "EvaluatorWrapper.h"
//#include "ChargeEvaluatorEwald.h"
using namespace std;
namespace py = boost::python;
const string RLCutType = "RLCut";



FixRLCut::FixRLCut(boost::shared_ptr<State> state_, string handle_)
    : FixPair(state_, handle_, "all", RLCutType, true, false, 1),
    epsHandle("eps"), sigHandle("sig"), rCutHandle("rCut"), alphaHandle("alpha"), r_aHandle("r_a"), a_gHandle("a_g"), sig_gHandle("sig_g"), r_gHandle("r_g")
{

    initializeParameters(epsHandle, epsilons);
    initializeParameters(sigHandle, sigmas);
    initializeParameters(rCutHandle, rCuts);
    initializeParameters(alphaHandle, alphas);
    initializeParameters(r_aHandle, r_as);
    initializeParameters(a_gHandle, a_gs);
    initializeParameters(sig_gHandle, sig_gs);
    initializeParameters(r_gHandle, r_gs);
    paramOrder = {rCutHandle, epsHandle, sigHandle, alphaHandle, r_aHandle, a_gHandle, sig_gHandle, r_gHandle};
    readFromRestart();
    canAcceptChargePairCalc = true;
    setEvalWrapper();
}

void FixRLCut::compute(int virialMode) {
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

void FixRLCut::singlePointEng(float *perParticleEng) {
    int nAtoms = state->atoms.size();
    int numTypes = state->atomParams.numTypes;
    GPUData &gpd = state->gpd;
    GridGPU &grid = state->gridGPU;
    int activeIdx = gpd.activeIdx();
    uint16_t *neighborCounts = grid.perAtomArray.d_data.data();
    float *neighborCoefs = state->specialNeighborCoefs;
    evalWrap->energy(nAtoms, gpd.xs(activeIdx), perParticleEng, neighborCounts, grid.neighborlist.data(), grid.perBlockArray.d_data.data(), state->devManager.prop.warpSize, paramsCoalesced.data(), numTypes, state->boundsGPU, neighborCoefs[0], neighborCoefs[1], neighborCoefs[2], gpd.qs(activeIdx), chargeRCut);

}

void FixRLCut::setEvalWrapper() {
    if (evalWrapperMode == "offload") {
        EvaluatorRL eval;
        evalWrap = pickEvaluator<EvaluatorRL, 8, true>(eval, chargeCalcFix);
    } else if (evalWrapperMode == "self") {
        EvaluatorRL eval;
        evalWrap = pickEvaluator<EvaluatorRL, 8, true>(eval, nullptr);
    }
}

bool FixRLCut::prepareForRun() {
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
    auto processAlpha = [] (float a) {
        return 1.0f/a;
    };
    auto processR_a = [] (float a) {
        return a;
    };
    auto processA_g = [] (float a) {
        return a;
    };
    auto processSig_g = [] (float a) {
        return 1.0f / (a * a);
    };
    auto processR_g = [] (float a) {
        return a;
    };
    prepareParameters(epsHandle, fillEps, processEps, false);
    prepareParameters(sigHandle, fillSig, processSig, false);
    prepareParameters(rCutHandle, fillRCut, processRCut, true, fillRCutDiag);
    prepareParameters(alphaHandle, fillEps, processAlpha, false);
    prepareParameters(r_aHandle, fillEps, processR_a, false);
    prepareParameters(a_gHandle, fillEps, processA_g, false);
    prepareParameters(sig_gHandle, fillEps, processSig_g, false);
    prepareParameters(r_gHandle, fillEps, processR_g, false);

    sendAllToDevice();
    setEvalWrapper();
    return true;
}

string FixRLCut::restartChunk(string format) {
    stringstream ss;
    ss << restartChunkPairParams(format);
    return ss.str();
}


bool FixRLCut::postRun() {

    return true;
}

void FixRLCut::addSpecies(string handle) {
    initializeParameters(epsHandle, epsilons);
    initializeParameters(sigHandle, sigmas);
    initializeParameters(rCutHandle, rCuts);
    initializeParameters(alphaHandle, alphas);
    initializeParameters(r_aHandle, r_as);
    initializeParameters(a_gHandle, a_gs);
    initializeParameters(sig_gHandle, sig_gs);
    initializeParameters(r_gHandle, r_gs);

}

vector<float> FixRLCut::getRCuts() {
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

void export_FixRLCut() {
    py::class_<FixRLCut, boost::shared_ptr<FixRLCut>, py::bases<FixPair>, boost::noncopyable > (
        "FixRLCut",
        py::init<boost::shared_ptr<State>, string> (py::args("state", "handle"))
    )
      ;
}
