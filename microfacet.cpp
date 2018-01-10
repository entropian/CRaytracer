#include "util/math.h"
#include "microfacet.h"

float MicrofacetDistribution_D(const vec3 wh, const MicrofacetDistribution* distrib)
{
    if(distrib->type == BECKMANN)
    {
        /*
          Float tan2Theta = Tan2Theta(wh);
          if (std::isinf(tan2Theta)) return 0.;
          Float cos4Theta = Cos2Theta(wh) * Cos2Theta(wh);
          return std::exp(-tan2Theta * (Cos2Phi(wh) / (alphax * alphax) +
          Sin2Phi(wh) / (alphay * alphay))) /
          (Pi * alphax * alphay * cos4Theta);
         */
        float tan_2_theta = tan2Theta(wh);
        if(isinf(tan_2_theta)) return 0.0f;
        float cos_4_theta = cos2Theta(wh) * cos2Theta(wh);
        return exp(-tan_2_theta * (cos2Phi(wh) / (distrib->alphax * distrib->alphax) +
                                   sin2Phi(wh) / (distrib->alphay * distrib->alphay))) /
            (PI * distrib->alphax * distrib->alphay * cos_4_theta);
    }else if(distrib->type == TROWBRIDGEREITZ)
    {
        float tan_2_theta = tan2Theta(wh);
        if(isinf(tan_2_theta)) return 0.0f;
        float cos_4_theta = cos2Theta(wh) * cos2Theta(wh);
        float e = (cos2Phi(wh) / (distrib->alphax * distrib->alphay) +
                   sin2Phi(wh) / (distrib->alphax * distrib->alphay)) * tan_2_theta;
        return 1.0f / (PI * distrib->alphax * distrib->alphay * cos_4_theta * (1.0f + e) * (1.0f + e));
    }
}

float MicrofacetDistribution_Lambda(const vec3 w, const MicrofacetDistribution* distrib)
{
    if(distrib->type == BECKMANN)
    {
        /*
          Float absTanTheta = std::abs(TanTheta(w));
          if (std::isinf(absTanTheta)) return 0.;
          // Compute _alpha_ for direction _w_
          Float alpha =
          std::sqrt(Cos2Phi(w) * alphax * alphax + Sin2Phi(w) * alphay * alphay);
          Float a = 1 / (alpha * absTanTheta);
          if (a >= 1.6f) return 0;
          return (1 - 1.259f * a + 0.396f * a * a) / (3.535f * a + 2.181f * a * a);
         */
        float abs_tan_theta = fabs(tanTheta(w));
        if(isinf(abs_tan_theta)) return 0.0f;
        float alpha = sqrtf(cos2Phi(w) * distrib->alphax * distrib->alphax +
                           sin2Phi(w) * distrib->alphay * distrib->alphay);
        float a = 1.0f / (alpha * abs_tan_theta);
        if(a >= 1.6f)
            return 0.0f;
        return (1.0f - 1.259f * a + 0.396f * a * a) /
            (3.535f * a + 2.181f * a * a);
    }else if(distrib->type == TROWBRIDGEREITZ)
    {
        float abs_tan_theta = fabs(tanTheta(w));
        if(isinf(abs_tan_theta)) return 0.0f;
        float alpha = sqrt(cos2Phi(w) * distrib->alphax * distrib->alphax +
                           sin2Phi(w) * distrib->alphay * distrib->alphay);
        float alpha_2_tan_2_theta = (alpha * abs_tan_theta) * (alpha * abs_tan_theta);
        return (-1.0f + sqrt(1.0f + alpha_2_tan_2_theta)) / 2.0f;
    }
}

float MicrofacetDistribution_G1(const vec3 w, const MicrofacetDistribution* distrib)
{
    return 1.0f / (1.0f + MicrofacetDistribution_Lambda(w, distrib));
}

float MicrofacetDistribution_G(const vec3 wo, const vec3 wi, const MicrofacetDistribution* distrib)
{
    return 1.0f / (1.0f + MicrofacetDistribution_Lambda(wo, distrib) + MicrofacetDistribution_Lambda(wi, distrib));
}

void MicrofacetDistribution_sample_wh(vec3 wh, const vec3 wo, const vec2 sample,
                                      const MicrofacetDistribution* distrib)
{
    /*
        // Compute $\tan^2 \theta$ and $\phi$ for Beckmann distribution sample
        Float tan2Theta, phi;
        if (alphax == alphay) {
            Float logSample = std::log(u[0]);
            if (std::isinf(logSample)) logSample = 0;
            tan2Theta = -alphax * alphax * logSample;
            phi = u[1] * 2 * Pi;
        } else {
            // Compute _tan2Theta_ and _phi_ for anisotropic Beckmann
            // distribution
            Float logSample = std::log(u[0]);
            if (std::isinf(logSample)) logSample = 0;
            phi = std::atan(alphay / alphax *
                            std::tan(2 * Pi * u[1] + 0.5f * Pi));
            if (u[1] > 0.5f) phi += Pi;
            Float sinPhi = std::sin(phi), cosPhi = std::cos(phi);
            Float alphax2 = alphax * alphax, alphay2 = alphay * alphay;
            tan2Theta = -logSample /
                        (cosPhi * cosPhi / alphax2 + sinPhi * sinPhi / alphay2);
        }
        // Map sampled Beckmann angles to normal direction _wh_
        Float cosTheta = 1 / std::sqrt(1 + tan2Theta);
        Float sinTheta = std::sqrt(std::max((Float)0, 1 - cosTheta * cosTheta));
        Vector3f wh = SphericalDirection(sinTheta, cosTheta, phi);
        if (!SameHemisphere(wo, wh)) wh = -wh;
        return wh;
     */
    float tan_2_theta, phi;
    if(distrib->alphax == distrib->alphay)
    {
        float log_sample = logf(sample[0]);
        if(isinf(log_sample))
            log_sample = 0.0f;
        tan_2_theta = -distrib->alphax * distrib->alphax * log_sample;
        phi = sample[1] * 2.0f * PI;
    }else
    {
        float log_sample = logf(sample[0]);
        if(isinf(log_sample))
            log_sample = 0.0f;
        phi = atanf(distrib->alphay / distrib->alphax *
                    tanf(2.0f * PI * sample[1] + 0.5f * PI));
        if(sample[1] > 0.5f) phi += PI;
        float sin_phi = sinf(phi), cos_phi = cosf(phi);
        float alphax2 = distrib->alphax * distrib->alphax;
        float alphay2 = distrib->alphay * distrib->alphay;
        tan_2_theta = -log_sample /
            (cos_phi * cos_phi / alphax2 + sin_phi * sin_phi / alphay2);
    }
    float cos_theta = 1.0f / sqrtf(1.0f + tan_2_theta);
    float sin_theta = sqrtf(max(0.0f, 1.0f - cos_theta * cos_theta));
    sphericalDirection(wh, sin_theta, cos_theta, phi);
    if(!sameHemisphere(wo, wh))
        vec3_negate(wh, wh);
}

float MicrofacetDistribution_pdf(const vec3 wo, const vec3 wh, const MicrofacetDistribution* distrib)
{
    return MicrofacetDistribution_D(wh, distrib) * absCosTheta(wh);
}
