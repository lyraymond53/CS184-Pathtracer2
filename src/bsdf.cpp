#include "bsdf.h"

#include <iostream>
#include <algorithm>
#include <utility>

using std::min;
using std::max;
using std::swap;

namespace CGL {

void make_coord_space(Matrix3x3& o2w, const Vector3D& n) {
  Vector3D z = Vector3D(n.x, n.y, n.z);
  Vector3D h = z;
  if (fabs(h.x) <= fabs(h.y) && fabs(h.x) <= fabs(h.z)) h.x = 1.0;
  else if (fabs(h.y) <= fabs(h.x) && fabs(h.y) <= fabs(h.z)) h.y = 1.0;
  else h.z = 1.0;

  z.normalize();
  Vector3D y = cross(h, z);
  y.normalize();
  Vector3D x = cross(z, y);
  x.normalize();

  o2w[0] = x;
  o2w[1] = y;
  o2w[2] = z;
}

// Mirror BSDF //

Spectrum MirrorBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum MirrorBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO: 1.2
  // Using BSDF::reflect(), implement sample_f for a mirror surface

  *pdf = (float) 1;

  reflect(wo, wi);

  return reflectance / (float) abs_cos_theta(*wi);
}

// Microfacet BSDF //

double MicrofacetBSDF::G(const Vector3D& wo, const Vector3D& wi) {
    return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
}

double MicrofacetBSDF::D(const Vector3D& h) {
  // TODO: 2.2
  // Compute Beckmann normal distribution function (NDF) here.
  // You will need the roughness alpha.

  double theta = getTheta(h);
  double alpha2 = alpha * alpha;

  double e = - ((tan(theta)*tan(theta))/alpha2);

  double den = PI * alpha2 * pow(cos_theta(h), 4);

  return exp(e)/den;
}

Spectrum MicrofacetBSDF::F(const Vector3D& wi) {
  // TODO: 2.3
  // Compute Fresnel term for reflection on dielectric-conductor interface.
  // You will need both eta and etaK, both of which are Spectrum.

  Spectrum eta_sq_k_sq = (eta * eta) + (k * k);
  Spectrum two_eta_cos = 2*eta*cos_theta(wi);
  double cos_sq = cos_theta(wi)*cos_theta(wi);

  Spectrum Rs = (eta_sq_k_sq - two_eta_cos + cos_sq) / (eta_sq_k_sq + two_eta_cos + cos_sq);
  Spectrum Rp = ((eta_sq_k_sq * cos_sq) - two_eta_cos + 1) / ((eta_sq_k_sq * cos_sq) + two_eta_cos + 1);

  Spectrum F = (Rs + Rp) / 2;

  return F;
}

Spectrum MicrofacetBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  // TODO: 2.1
  // Implement microfacet model here

  if ((wo.z <= 0) || (wi.z <= 0)) {
      return Spectrum();
  }

  Vector3D h = (wo + wi).unit();

  Spectrum a = (F(wi) * G(wo, wi) * D(h)) / (4 * wo.z * wi.z);

  return a;
}

Spectrum MicrofacetBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO: 2.4
  // *Importance* sample Beckmann normal distribution function (NDF) here.
  // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
  //       and return the sampled BRDF value.

//  *wi = cosineHemisphereSampler.get_sample(pdf); //placeholder

  Vector2D r = sampler.get_sample();
  float alpha2 = alpha*alpha;

  float thetah = (float) atan(sqrt(-alpha2 * log(1.0 - r.x)));
    float tanthetah = tan(thetah);
    float sinthetah = sin(thetah);
    float costhetah = cos(thetah);


  float phih = 2.0 * PI * r.y;
    float cosphi = cos(phih);
    float sinphi = sin(phih);

  float exponent = exp(-(tanthetah * tanthetah)/alpha2);

  float ptheta = (float) (2.0 * sinthetah * exponent) / (float) (alpha2 * pow(costhetah, 3));
  float pphi = (float) 1.0 / (float) (2.0 * PI);

  Vector3D h = Vector3D((sinthetah * cosphi), (sinphi * sinthetah), (costhetah));

  *wi = (2 * dot(wo, h) * h) - wo;

  float pwh = (ptheta * pphi) / sinthetah;

  if (wi->z <= 0) {
      *pdf = 0;
      return Spectrum();
  }

  *pdf = (pwh / (float) (4.0 * dot(*wi, h)));

  return MicrofacetBSDF::f(wo, *wi);
}

// Refraction BSDF //

Spectrum RefractionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum RefractionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  return Spectrum();
}

// Glass BSDF //

Spectrum GlassBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum GlassBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO: 1.4
  // Compute Fresnel coefficient and either reflect or refract based on it.
  float eta;

  if (refract(wo, wi, ior) == false) {
      reflect(wo, wi);
      *pdf = 1;
      return reflectance / (float) abs_cos_theta(*wi);
  }
  else {
      float r0 = (1 - ior)/(1 + ior);
      float schlick = r0*r0 + (1-r0*r0)*pow((1 - abs_cos_theta(*wi)), 5);

      if (coin_flip(schlick)) {
          reflect(wo, wi);
          *pdf = schlick;
          return schlick * reflectance / (float) abs_cos_theta(*wi);
      }
      else {
          refract(wo, wi, ior);
          *pdf = 1 - schlick;
          if (wo.z > 0) {
              eta = 1/ior;
          } else {
              eta = ior;
          }
          return (1 - schlick) * transmittance / (float) abs_cos_theta(*wi) / (eta*eta);
      }
  }
}

void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {

  // TODO: 1.1
  // Implement reflection of wo about normal (0,0,1) and store result in wi.

  *wi = Vector3D(-wo.x, -wo.y, wo.z);
}

bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior) {

  // TODO: 1.3
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
  // ray entering the surface through vacuum.
  float n;
  float temp;

  if (wo.z > 0) {
      n = 1/ior;
      temp = -1;
  } else {
      n = ior;
      temp = 1;
  }

  float internal = 1 - n*n * (1 - wo.z*wo.z);

  if (internal < 0) {
      return false;
  }
  *wi = Vector3D(-n*wo.x, -n*wo.y, temp * sqrt(1-n*n*(1-wo.z*wo.z)));
  return true;
}

// Emission BSDF //

Spectrum EmissionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum EmissionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *pdf = 1.0 / PI;
  *wi  = sampler.get_sample(pdf);
  return Spectrum();
}

} // namespace CGL
