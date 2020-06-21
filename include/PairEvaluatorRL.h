#pragma once
#ifndef EVALUATOR_RL
#define EVALUATOR_RL

#include "cutils_math.h" 

class EvaluatorRL {
	public:
		inline __device__ float3 force(float3 dr, float params[8], float lenSqr, float multiplier) {
			if (multiplier) {
				float eps = params[1];
				float sig64 = params[2];
				float alpha_inv = params[3];
				float r_a = params[4];
				float a_g = params[5];
				float sig_g_sqr_inv = params[6];
				float r_g = params[7];

				float r2inv = 1.0f / lenSqr;
				float r1inv = sqrt(r2inv);
				float r4inv = r2inv * r2inv; 
				float p1 = 64.0f * sig64 * r4inv * r4inv * r4inv * r1inv;
				float p2 = alpha_inv * exp(alpha_inv * -1 * ((1.0f/r1inv) - r_a));
				float p3 = a_g * sig_g_sqr_inv * (1.0f/r1inv - r_g) * exp(-1.0f * 2.0f * sig_g_sqr_inv * (1.0f/r1inv - r_g) * (1.0f/r1inv - r_g));
				
				float forceScalar = (p1 - p2 + p3) * eps * multiplier;
				return dr * forceScalar;
			}
			return make_float3(0, 0, 0);
		}
		inline __device__ float energy(float params[8], float lenSqr, float multiplier) {
			if (multiplier) {
				float eps = params[1];
				float sig64 = params[2];
				float alpha_inv = params[3];
				float r_a = params[4];
				float a_g = params[5];
				float sig_g_sqr_inv = params[6];
				float r_g = params[7];

				float r2inv = 1.0f / lenSqr;
				float r1inv = sqrt(r2inv);
				float r4inv = r2inv * r2inv; 
				float p1 = sig64 * r4inv * r4inv * r4inv;
				float p2 = exp(-1.0f * alpha_inv * ((1.0f/r1inv) - r_a));
				float p3 = a_g * exp(-.5f * sig_g_sqr_inv * (1.0f/r1inv - r_g) * (1.0f/r1inv - r_g));

				float rCutSqr = params[0];
				float rCut6 = rCutSqr * rCutSqr * rCutSqr;
				float sig4 = std::pow(sig64, 1/3);
				float sig2 = sqrt(sig4);
				float sig6InvRCut6 = sig4 * sig2 / rCut6;
				float offsetOver4Eps = sig6InvRCut6 * (sig6InvRCut6 - 1.0f);

				return 0.5f * eps * (p1 - p2 + p3 - (offsetOver4Eps * 4)) * multiplier ;
			}
			return 0;
		}
};

#endif
