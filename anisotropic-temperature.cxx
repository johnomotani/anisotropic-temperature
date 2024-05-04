#include <bout/physicsmodel.hxx>

// This macro should already have been defined in bout/field.hxx, but it seems not to be
// useable in this file, so copy-paste the definition in here.
/*!
 * This macro takes a function \p func, which is
 * assumed to operate on a single BoutReal and return
 * a single BoutReal, and wraps it up into a function
 * of a Field called \p name.
 *
 * @param name  The name of the function to define
 * @param func  The function to apply to each value
 *
 * If CHECK >= 1, checks if the Field is allocated
 *
 * Loops over the entire domain, applies function,
 * and uses checkData() to, if CHECK >= 3, check
 * result for non-finite numbers
 *
 */
#ifdef FIELD_FUNC
#error This macro has already been defined
#else
#define FIELD_FUNC(name, func)                                     \
  template <typename T, typename = bout::utils::EnableIfField<T>>  \
  inline T name(const T& f, const std::string& rgn = "RGN_ALL") {  \
    AUTO_TRACE();                                                  \
    /* Check if the input is allocated */                          \
    checkData(f);                                                  \
    /* Define and allocate the output result */                    \
    T result{emptyFrom(f)};                                        \
    BOUT_FOR(d, result.getRegion(rgn)) { result[d] = func(f[d]); } \
    checkData(result);                                             \
    return result;                                                 \
  }
#endif

class AnisotropicTemperature : public PhysicsModel {

  int init(bool UNUSED(restarting)) {

    solver->add(n, "n");
    solver->add(V_i, "V_i");
    solver->add(ppar_i, "ppar_i");
    solver->add(pperp_i, "pperp_i");

    solver->add(ppar_e, "ppar_e");
    solver->add(pperp_e, "pperp_e");

    return 0;
  }

  int rhs(BoutReal UNUSED(t)) {

    mesh->communicate(n, V_i, ppar_i, pperp_i, ppar_e, pperp_e);

    Field3D Epar = - Grad_par(ppar_e) / n;

    Field3D Tperp_i = pperp_i / n;
    Field3D Tpar_i = ppar_i / n;
    Field3D Tperp_e = pperp_e / n;
    Field3D Tpar_e = ppar_e / n;

    Field3D nu_ee = n / Tperp_e / sqrt(Tpar_e);
    Field3D nu_ei = nu_ee;
    Field3D nu_ii = n / Tperp_i / sqrt(m_i * Tpar_i);
    Field3D nu_ie = nu_ee / m_i;

    Field3D betapar_i = m_i / Tpar_i;
    Field3D betapar_e = 1.0 / Tpar_e;
    Field3D betapar_ii = 0.5 * betapar_i;
    Field3D betapar_ie = betapar_i * betapar_e / (betapar_i + betapar_e);
    Field3D betapar_ee = 0.5 * betapar_e;
    Field3D betapar_ei = betapar_ie;
    Field3D betaperp_i = m_i / Tperp_i;
    Field3D betaperp_e = 1.0 / Tperp_e;
    Field3D betaperp_ii = 0.5 * betaperp_i;
    Field3D betaperp_ie = betaperp_i * betaperp_e / (betaperp_i + betaperp_e);
    Field3D betaperp_ee = 0.5 * betaperp_e;
    Field3D betaperp_ei = betaperp_ie;

    Field3D alpha_ii = betapar_ii / betaperp_ii;
    Field3D alpha_ie = betapar_ie / betaperp_ie;
    Field3D alpha_ee = betapar_ee / betaperp_ee;
    Field3D alpha_ei = alpha_ie;

    Field3D X_ii = alpha_ii - 1.0;
    Field3D X_ie = alpha_ie - 1.0;
    Field3D X_ee = alpha_ee - 1.0;
    Field3D X_ei = X_ie;

    Field3D Kii_004 = K004(X_ii);
    Field3D Kii_202 = K202(X_ii);
    Field3D Kii_200 = K200(X_ii);
    Field3D Kii_002 = K002(X_ii);
    Field3D Kii_220 = K220(X_ii);

    Field3D Kie_200 = K200(X_ie);
    Field3D Kie_002 = K002(X_ie);

    Field3D Kee_004 = K004(X_ee);
    Field3D Kee_202 = K202(X_ee);
    Field3D Kee_200 = K200(X_ee);
    Field3D Kee_002 = K002(X_ee);
    Field3D Kee_220 = K220(X_ee);
    Field3D Kee_222 = K222(X_ee);
    Field3D Kee_204 = K204(X_ee);
    Field3D Kee_006 = K006(X_ee);

    Field3D Kei_200 = Kie_200;
    Field3D Kei_002 = Kie_002;

    Field3D psi_i = SQ(alpha_ii) * (Kii_004 - Kii_202) + 0.5 * alpha_ii * (Kii_200 - Kii_002);
    Field3D phi_i = 4.0 * Kii_220 - 2.0 * Kii_202 - Kii_200 + Kii_002;
    Field3D psi_e = SQ(alpha_ee) * (Kee_004 - Kee_202) + 0.5 * alpha_ee * (Kee_200 - Kee_002);
    Field3D phi_e = 4.0 * Kee_220 - 2.0 * Kee_202 - Kee_200 + Kee_002;

    Field3D cperp_i = - 4.0 * nu_ii * (6.0 * alpha_ii * Kii_202 + 0.5 * phi_i)
                      - 4.0 * nu_ie * (2.0 * Kie_200 + alpha_ie * Kie_002);
    Field3D cpar_i = 2.0 * psi_i * nu_ii;
    Field3D eperp_i = 12.0 * phi_i * nu_ii;
    Field3D epar_i = -12.0 * psi_i * nu_ii - 12.0 * nu_ie * alpha_ie * Kie_002;

    Field3D cperp_e = -4.0 * nu_ee * (6.0 * alpha_ee * Kee_202 + 0.5 * phi_e)
                      + 4.0 * alpha_ee * nu_ei * (- 32.0 * Kee_222 + 4.0 * Kee_204
                                                  + 10.0 * Kee_202 - 2.0 * Kee_004
                                                  - Kee_002);
    Field3D cpar_e = 2.0 * nu_ee * psi_e
                     + 4.0 * SQ(alpha_ee) * nu_ei
                       * (-8.0/3.0 * alpha_ee * Kee_204 + alpha_ee * Kee_006
                          + 4.0 * Kee_202 - (1 - alpha_ee / 3.0) * Kee_004
                          - 0.5 * Kee_002);
    Field3D eperp_e = 12.0 * nu_ee * phi_e
                      + 12.0 * nu_ei * (16.0 * alpha_ee * Kee_222 - 4.0 * alpha_ee * Kee_204
                                        - 2.0 * (2.0 * alpha_ee - 1.0) * Kee_202
                                        + 2.0 * alpha_ee * Kee_004 - Kee_002);
    Field3D epar_e = - 12.0 * nu_ee * psi_e
                     + 4.0 * nu_ei * alpha_ee * (4.0 * SQ(alpha_ee) * Kee_204
                                                 - 2.0 * SQ(alpha_ee) * Kee_006
                                                 - 6.0 * alpha_ee * Kee_202
                                                 + 4.0 * alpha_ee * Kee_004
                                                 - 1.5 * Kee_002);

    Field3D grad_Tperp_i = Grad_par(pperp_i/n);
    Field3D grad_Tpar_i = Grad_par(ppar_i/n);
    Field3D grad_Tperp_e = Grad_par(pperp_e/n);
    Field3D grad_Tpar_e = Grad_par(ppar_e/n);

    Field3D S_perp_i_par = S_perp_s_par(cperp_i, cpar_i, eperp_i, epar_i, ppar_i,
                                        grad_Tperp_i, grad_Tpar_i, m_i);
    Field3D S_par_i_par = S_par_s_par(cperp_i, cpar_i, eperp_i, epar_i, ppar_i,
                                      grad_Tperp_i, grad_Tpar_i, m_i);
    Field3D S_perp_e_par = S_perp_s_par(cperp_e, cpar_e, eperp_e, epar_e, ppar_e,
                                        grad_Tperp_e, grad_Tpar_e, 1.0);
    Field3D S_par_e_par = S_par_s_par(cperp_e, cpar_e, eperp_e, epar_e, ppar_e,
                                      grad_Tperp_e, grad_Tpar_e, 1.0);

    Field3D J_ii_perp_perp = J_sr_perp_perp(n, nu_ii, m_i, m_i, Kii_200, Tperp_i, Tperp_i,
                                            Kii_002);
    Field3D J_ie_perp_perp = J_sr_perp_perp(n, nu_ie, m_i, 1.0, Kie_200, Tperp_i, Tperp_e,
                                            Kie_002);
    Field3D J_ii_par_par = J_sr_par_par(n, nu_ii, m_i, m_i, alpha_ii, Kii_002, Tpar_i,
                                        Tpar_i, Kii_200);
    Field3D J_ie_par_par = J_sr_par_par(n, nu_ie, m_i, 1.0, alpha_ie, Kie_002, Tpar_i,
                                        Tpar_e, Kie_200);
    Field3D J_ee_perp_perp = J_sr_perp_perp(n, nu_ee, 1.0, 1.0, Kee_200, Tperp_e, Tperp_e,
                                            Kee_002);
    Field3D J_ei_perp_perp = J_sr_perp_perp(n, nu_ei, 1.0, m_i, Kei_200, Tperp_e, Tperp_i,
                                            Kei_002);
    Field3D J_ee_par_par = J_sr_par_par(n, nu_ee, 1.0, 1.0, alpha_ee, Kee_002, Tpar_e,
                                        Tpar_e, Kee_200);
    Field3D J_ei_par_par = J_sr_par_par(n, nu_ei, 1.0, m_i, alpha_ei, Kei_002, Tpar_e,
                                        Tpar_i, Kei_200);


    ddt(n) = - Div_par_flux(V_i, n);

    ddt(V_i) = - Vpar_Grad_par(V_i, V_i)
               - Grad_par(ppar_i) / (m_i * n)
               - Epar / m_i;

    //ddt(pperp_i) = - Vpar_Grad_par(V_i, pperp_i)
    //               - pperp_i * Grad_par(V_i)
    //               - Grad_par(S_perp_i_par)
    //               + J_ii_perp_perp + J_ie_perp_perp;
    ddt(pperp_i) = 0.0;

    //ddt(ppar_i) = - Vpar_Grad_par(V_i, ppar_i)
    //              - 3.0 * ppar_i * Grad_par(V_i)
    //              - Grad_par(S_par_i_par)
    //              + J_ii_par_par + J_ie_par_par;
    ddt(ppar_i) = 0.0;

    //ddt(pperp_e) = - Vpar_Grad_par(V_i, pperp_e)
    //               - pperp_e * Grad_par(V_i)
    //               - Grad_par(S_perp_e_par)
    //               + J_ee_perp_perp + J_ei_perp_perp;
    ddt(pperp_e) = 0.0;

    ddt(ppar_e) = - Vpar_Grad_par(V_i, ppar_e)
                  - 3.0 * ppar_e * Grad_par(V_i)
                  //- Grad_par(S_par_e_par)
                  //+ J_ee_par_par + J_ei_par_par
                  ;

    return 0;
  }

  Field3D n;
  Field3D V_i, ppar_i, pperp_i;
  Field3D ppar_e, pperp_e;

  // Deuteron mass normalised to electron mass
  const BoutReal m_i = 3.3435837724e-27 / 9.1093837015e-31;

  const BoutReal epsilon_near_zero = 1.0e-6;

  BoutReal phi(const BoutReal& X) {
    // copysign means the result returned has the same sign as X. Note that this function
    // should never be called very near to X=0.
    return copysign(std::atan(std::sqrt(std::abs(X))) / std::sqrt(std::abs(X)), X);
  }

  BoutReal K004_singular(const BoutReal& X) {
    return (2.0 + 1.0 / (1.0 + X) - 3.0 * phi(X)) / SQ(X);
  }
  BoutReal K004(const BoutReal& X) {
    // We know that the function is non-singular and linear in X near X=0. To avoid
    // numerical issues, replace with an approximate version which is continuous at
    // +/-epsilon_near_zero when X is near 0.
    if (abs(X) < epsilon_near_zero) {
      return (K004_singular(-epsilon_near_zero) * (epsilon_near_zero - X)
              + K004_singular(epsilon_near_zero) * (X + epsilon_near_zero))
             / (2.0 * epsilon_near_zero);
    } else {
      return K004_singular(X);
    }
  }
  FIELD_FUNC(K004, K004)

  BoutReal K202_singular(const BoutReal& X) {
    return 0.5 * (-3.0 + (3.0 * X) * phi(X)) / SQ(X);
  }
  BoutReal K202(const BoutReal& X) {
    // We know that the function is non-singular and linear in X near X=0. To avoid
    // numerical issues, replace with an approximate version which is continuous at
    // +/-epsilon_near_zero when X is near 0.
    if (abs(X) < epsilon_near_zero) {
      return (K202_singular(-epsilon_near_zero) * (epsilon_near_zero - X)
              + K202_singular(epsilon_near_zero) * (X + epsilon_near_zero))
             / (2.0 * epsilon_near_zero);
    } else {
      return K202_singular(X);
    }
  }
  FIELD_FUNC(K202, K202)

  BoutReal K200_singular(const BoutReal& X) {
    return (-1.0 + (1.0 + X) * phi(X)) / X;
  }
  BoutReal K200(const BoutReal& X) {
    // We know that the function is non-singular and linear in X near X=0. To avoid
    // numerical issues, replace with an approximate version which is continuous at
    // +/-epsilon_near_zero when X is near 0.
    if (abs(X) < epsilon_near_zero) {
      return (K200_singular(-epsilon_near_zero) * (epsilon_near_zero - X)
              + K200_singular(epsilon_near_zero) * (X + epsilon_near_zero))
             / (2.0 * epsilon_near_zero);
    } else {
      return K200_singular(X);
    }
  }
  FIELD_FUNC(K200, K200)

  BoutReal K002_singular(const BoutReal& X) {
    return 2.0 * (1.0 - phi(X)) / X;
  }
  BoutReal K002(const BoutReal& X) {
    // We know that the function is non-singular and linear in X near X=0. To avoid
    // numerical issues, replace with an approximate version which is continuous at
    // +/-epsilon_near_zero when X is near 0.
    if (abs(X) < epsilon_near_zero) {
      return (K002_singular(-epsilon_near_zero) * (epsilon_near_zero - X)
              + K002_singular(epsilon_near_zero) * (X + epsilon_near_zero))
             / (2.0 * epsilon_near_zero);
    } else {
      return K002_singular(X);
    }
  }
  FIELD_FUNC(K002, K002)

  BoutReal K220_singular(const BoutReal& X) {
    return 0.125 * (3.0 + X + (1.0 + X) * (X - 3.0) * phi(X)) / SQ(X);
  }
  BoutReal K220(const BoutReal& X) {
    // We know that the function is non-singular and linear in X near X=0. To avoid
    // numerical issues, replace with an approximate version which is continuous at
    // +/-epsilon_near_zero when X is near 0.
    if (abs(X) < epsilon_near_zero) {
      return (K220_singular(-epsilon_near_zero) * (epsilon_near_zero - X)
              + K220_singular(epsilon_near_zero) * (X + epsilon_near_zero))
             / (2.0 * epsilon_near_zero);
    } else {
      return K220_singular(X);
    }
  }
  FIELD_FUNC(K220, K220)

  BoutReal K222_singular(const BoutReal& X) {
    return 0.0625 * (15.0 + X + (-15.0 - 6.0 * X + SQ(X)) * phi(X)) / (X*X*X);
  }
  BoutReal K222(const BoutReal& X) {
    // We know that the function is non-singular and linear in X near X=0. To avoid
    // numerical issues, replace with an approximate version which is continuous at
    // +/-epsilon_near_zero when X is near 0.
    if (abs(X) < epsilon_near_zero) {
      return (K222_singular(-epsilon_near_zero) * (epsilon_near_zero - X)
              + K222_singular(epsilon_near_zero) * (X + epsilon_near_zero))
             / (2.0 * epsilon_near_zero);
    } else {
      return K222_singular(X);
    }
  }
  FIELD_FUNC(K222, K222)

  BoutReal K204_singular(const BoutReal& X) {
    return 0.25 * (-13.0 - 2.0 / (1.0 + X) + (15.0 + 3.0 * X) * phi(X)) / (X*X*X);
  }
  BoutReal K204(const BoutReal& X) {
    // We know that the function is non-singular and linear in X near X=0. To avoid
    // numerical issues, replace with an approximate version which is continuous at
    // +/-epsilon_near_zero when X is near 0.
    if (abs(X) < epsilon_near_zero) {
      return (K204_singular(-epsilon_near_zero) * (epsilon_near_zero - X)
              + K204_singular(epsilon_near_zero) * (X + epsilon_near_zero))
             / (2.0 * epsilon_near_zero);
    } else {
      return K204_singular(X);
    }
  }
  FIELD_FUNC(K204, K204)

  BoutReal K006_singular(const BoutReal& X) {
    return 0.5 * (8.0 + 9.0 / (1.0 + X) - 2.0 / SQ(1.0 + X) - 15.0 * phi(X)) / (X*X*X);
  }
  BoutReal K006(const BoutReal& X) {
    // We know that the function is non-singular and linear in X near X=0. To avoid
    // numerical issues, replace with an approximate version which is continuous at
    // +/-epsilon_near_zero when X is near 0.
    if (abs(X) < epsilon_near_zero) {
      return (K006_singular(-epsilon_near_zero) * (epsilon_near_zero - X)
              + K006_singular(epsilon_near_zero) * (X + epsilon_near_zero))
             / (2.0 * epsilon_near_zero);
    } else {
      return K006_singular(X);
    }
  }
  FIELD_FUNC(K006, K006)

  Field3D S_perp_s_par(const Field3D& cperp_s, const Field3D& cpar_s,
                       const Field3D& eperp_s, const Field3D& epar_s,
                       const Field3D& ppar_s, const Field3D& grad_Tperp_s,
                       const Field3D& grad_Tpar_s, const BoutReal& m_s) {
    return ppar_s / m_s * (epar_s * grad_Tperp_s - 3.0 * cpar_s * grad_Tpar_s)
                          / (cperp_s * epar_s - cpar_s * eperp_s);
  }

  Field3D S_par_s_par(const Field3D& cperp_s, const Field3D& cpar_s,
                      const Field3D& eperp_s, const Field3D& epar_s,
                      const Field3D& ppar_s, const Field3D& grad_Tperp_s,
                      const Field3D& grad_Tpar_s, const BoutReal& m_s) {
    return ppar_s / m_s * (eperp_s * grad_Tperp_s - 3.0 * cperp_s * grad_Tpar_s)
                          / (cpar_s * eperp_s - cperp_s * epar_s);
  }

  Field3D J_sr_perp_perp(const Field3D& n, const Field3D& nu_sr, const BoutReal& m_s,
                         const BoutReal& m_r, const Field3D& Ksr_200,
                         const Field3D& Tperp_s, const Field3D& Tperp_r,
                         const Field3D& Ksr_002) {
    Field3D betaperp_r = m_r / Tperp_r;
    Field3D betaperp_s = m_s / Tperp_s;
    Field3D betaperp_sr = betaperp_s * betaperp_r / (betaperp_s + betaperp_r);
    return 4.0 * m_s * n * nu_sr / (m_r + m_s) * (2.0 * Ksr_200 * (Tperp_r - Tperp_s)
                                                  + m_r / betaperp_sr * (Ksr_002 - Ksr_200));
  }

  Field3D J_sr_par_par(const Field3D& n, const Field3D& nu_sr, const BoutReal& m_s,
                       const BoutReal& m_r, const Field3D& alpha_sr,
                       const Field3D& Ksr_002, const Field3D& Tpar_s,
                       const Field3D& Tpar_r, const Field3D& Ksr_200) {
    Field3D betapar_r = m_r / Tpar_r;
    Field3D betapar_s = m_s / Tpar_s;
    Field3D betapar_sr = betapar_s * betapar_r / (betapar_s + betapar_r);
    return 8.0 * m_s * n * nu_sr / (m_r + m_s) * alpha_sr *
           (2.0 * Ksr_002 * (Tpar_r - Tpar_s) + m_r / betapar_sr * (Ksr_200 - Ksr_002));
  }
};

BOUTMAIN(AnisotropicTemperature);
