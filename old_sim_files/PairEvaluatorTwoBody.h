#pragma once
#ifndef EVALUATOR_TWOBODY
#define EVALUATOR_TWOBODY

#include "cutils_math.h" 

class EvaluatorTwoBody {
	public:
		inline __device__ float3 force(float3 dr, float params[9], float lenSqr, float multiplier) {
			if (multiplier) {
				float eps = params[1];
				float sig64 = params[2];
				float sigma_R = params[3];
				float r_A = params[4];
				float alpha_A_inv = params[5];
				float G = params[6];
				float r_G = params[7];
				float sigma_G_sqr_2 = params[8];

				float r1 = sqrt(lenSqr)
				float r2inv = 1.0f / lenSqr;
				float p1 = 64.0f * pow(sigma_R, 64) / pow(r1, 65); 
				float p2 = -1.0f * alpha_A_inv * exp((r_A-r1) * alpha_A_inv);
				float p3 = G * (r1 - r_G) * .5f * sigma_G_sqr_2 * exp(-1.0f * (r1 - r_G) * (r1 - r_G) * sigma_G_sqr_2)
				
				float forceScalar = (p1 - p2 + p3) * (eps / 24.0f) * multiplier;
				return dr * forceScalar;
			}
			return make_float3(0, 0, 0);
		}
		inline __device__ float energy(float params[9], float lenSqr, float multiplier) {
			if (multiplier) {
				float eps = params[1];
				float sig64 = params[2];
				float sigma_R = params[3];
				float r_A = params[4];
				float alpha_A_inv = params[5];
				float G = params[6];
				float r_G = params[7];
				float sigma_G_sqr_2 = params[8];

				float r1 = sqrt(lenSqr)
				float r2inv = 1.0f / lenSqr;
				float r4inv = r2inv * r2inv; 
				float p1 = pow((sigma_R / r1), 64);
				float p2 = exp(-1.0f * (r1 - r_A) * alpha_A_inv)
				float p3 = G * exp(-1.0f * (r1 - r_G) * (r1 - r_G) * sigma_G_sqr_2)

				float rCutSqr = params[0];
				float rCut6 = rCutSqr * rCutSqr * rCutSqr;
				float sig4 = std::pow(sig64, 1/3);
				float sig2 = sqrt(sig4);
				float sig6InvRCut6 = sig4 * sig2 / rCut6;
				float offsetOver4Eps = sig6InvRCut6 * (sig6InvRCut6 - 1.0f);

				return 0.5f * (eps / 24.0f) * (p1 - p2 + p3 - (offsetOver4Eps * 4)) * multiplier ;
			}
			return 0;
		}
};

#endif
