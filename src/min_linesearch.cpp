#include "min_linesearch.h"

using namespace KMC;
using namespace std;

// ALPHA_MAX = max alpha allowed to avoid long backtracks
// ALPHA_REDUCE = reduction ratio, should be in range [0.5,1)
// BACKTRACK_SLOPE, should be in range (0,0.5]
// QUADRATIC_TOL = tolerance on alpha0, should be in range [0.1,1)
// EMACH = machine accuracy limit of energy changes (1.0e-8)
// EPS_QUAD = tolerance for quadratic projection

const prec_t ALPHA_MAX = 1.0;
const prec_t ALPHA_REDUCE = 0.5;
const prec_t BACKTRACK_SLOPE = 0.4;
const prec_t QUADRATIC_TOL = 0.1;
//#define EMACH 1.0e-8
const prec_t EMACH = 1.0e-8;
const prec_t EPS_QUAD = 1.0e-28;

/* ---------------------------------------------------------------------- */

KMC::minLineSearch::minLineSearch():min(), g(nullptr), h(nullptr)
{
    searchflag = 1;
}

KMC::minLineSearch::minLineSearch(prec_t tol_ene, prec_t tol_force):min(tol_ene, tol_force), g(nullptr), h(nullptr)
{
    searchflag = 1;
}

/* ---------------------------------------------------------------------- */

KMC::minLineSearch::~minLineSearch()
{
    delete [] g;
    delete [] h;
}

/* ---------------------------------------------------------------------- */

void KMC::minLineSearch::init(potential *pot)
{
    min::init(pot);
    linemin = &minLineSearch::linemin_quadratic;
}

/* ---------------------------------------------------------------------- */

void KMC::minLineSearch::setup_style()
{
    nvec = system->get_atom_num() * 3;
    x0.resize(nvec);
    if (g)
        delete [] g;
    g = new prec_t [nvec];
    if (h)
        delete [] h;
    h = new prec_t [nvec];
}

/* ----------------------------------------------------------------------
   set current vector lengths and pointers
   called after atoms have migrated
------------------------------------------------------------------------- */

void KMC::minLineSearch::reset_vectors()
{
  // atomic dof

  if (system->get_atom_num()) xvec = system->positions[0];
  if (system->get_atom_num()) fvec = system->forces[0];
  if (system->get_atom_num()) mvec = const_cast<bool_t *>(system->get_moves()[0]);
}

/* ----------------------------------------------------------------------
   line minimization methods
   find minimum-energy starting at x along h direction
   input args:   eoriginal = energy at initial x
   input extra:  n,x,x0,f,h for atomic, extra global, extra per-atom dof
   output args:  return 0 if successful move, non-zero alpha
                 return non-zero if failed
                 alpha = distance moved along h for x at min eng config
                 update neval counter of eng/force function evaluations
                 output extra: if fail, energy_force() of original x
                 if succeed, energy_force() at x + alpha*h
                 atom->x = coords at new configuration
                 atom->f = force at new configuration
                 ecurrent = energy of new configuration
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   linemin: backtracking line search (Proc 3.1, p 41 in Nocedal and Wright)
   uses no gradient info, but should be very robust
   start at maxdist, backtrack until energy decrease is sufficient
------------------------------------------------------------------------- */

index_t KMC::minLineSearch::linemin_backtrack(prec_t orig_ene, prec_t &alpha)
{
    int i;
    prec_t fdothme, hme;
    prec_t de_ideal, de;

    // fdothall = projection of search dir along downhill gradient
    // if search direction is not downhill, exit with error

    fdothme = 0.0;
    for (i = 0; i < nvec; i++) fdothme += fvec[i] * h[i];
    if (fdothme <= 0.0) return DOWNHILL;

    // set alpha so no dof is changed by more than max allowed amount
    // for atom coords, max amount = dmax
    // for extra per-atom dof, max amount = extra_max[]
    // for extra global dof, max amount is set by fix
    // also ensure alpha <= ALPHA_MAX
    // else will have to backtrack from huge value when forces are tiny
    // if all search dir components are already 0.0, exit with error

    hme = 0.0;
    for (i = 0; i < nvec; i++) hme = max(hme,fabs(h[i]));
    alpha = std::min(ALPHA_MAX, dist_max/hme);
    if (hme == 0.0) return ZEROFORCE;

    // store box and values of all dof at start of linesearch

    // fix_minimize->store_box();
    for (i = 0; i < nvec; i++) x0[i] = xvec[i];

    // // important diagnostic: test the gradient against energy
    // prec_t etmp;
    // prec_t alphatmp = alpha*1.0e-4;
    // etmp = alpha_step(alphatmp,1);
    // printf("alpha = %g dele = %g dele_force = %g err = %g\n",
    //        alphatmp,etmp-eoriginal,-alphatmp*fdothall,
    //        etmp-eoriginal+alphatmp*fdothall);
    // alpha_step(0.0,1);

    // backtrack with alpha until energy decrease is sufficient

    while (true) {
        curr_ene = alpha_step(alpha);

        // if energy change is better than ideal, exit with success

        de_ideal = -BACKTRACK_SLOPE*alpha*fdothme;
        de = curr_ene - orig_ene;
        if (de <= de_ideal) 
            return 0;

        // reduce alpha

        alpha *= ALPHA_REDUCE;

        // backtracked too much
        // reset to starting point
        // if de is positive, exit with error
        // if de is negative, exit with ETOL

        if (alpha <= 0.0 || de_ideal >= -EMACH) {
            curr_ene = alpha_step(0.0);
            if (de < 0.0) return ETOL;
            else return ZEROALPHA;
        }
    }
}

/* ----------------------------------------------------------------------
    // linemin: quadratic line search (adapted from Dennis and Schnabel)
    // The objective function is approximated by a quadratic
    // function in alpha, for sufficiently small alpha.
    // This idea is the same as that used in the well-known secant
    // method. However, since the change in the objective function
    // (difference of two finite numbers) is not known as accurately
    // as the gradient (which is close to zero), all the expressions
    // are written in terms of gradients. In this way, we can converge
    // the LAMMPS forces much closer to zero.
    //
    // We know E, Eprev, fh, fhprev. The Taylor series about alpha_prev
    // truncated at the quadratic term is:
    //
    //     E = Eprev - del_alpha*fhprev + (1/2)del_alpha^2*Hprev
    //
    // and
    //
    //     fh = fhprev - del_alpha*Hprev
    //
    // where del_alpha = alpha-alpha_prev
    //
    // We solve these two equations for Hprev and E=Esolve, giving:
    //
    //     Esolve = Eprev - del_alpha*(f+fprev)/2
    //
    // We define relerr to be:
    //
    //      relerr = |(Esolve-E)/Eprev|
    //             = |1.0 - (0.5*del_alpha*(f+fprev)+E)/Eprev|
    //
    // If this is accurate to within a reasonable tolerance, then
    // we go ahead and use a secant step to fh = 0:
    //
    //      alpha0 = alpha - (alpha-alphaprev)*fh/delfh;
    //
------------------------------------------------------------------------- */

index_t KMC::minLineSearch::linemin_quadratic(prec_t orig_ene, prec_t &alpha)
{
    int i,m,n;
    prec_t fdothme,hme;
    prec_t de_ideal,de;
    prec_t delfh,engprev,relerr,alphaprev,fhprev,fh,alpha0;
    prec_t dot[2];
    prec_t alphamax;

    // fdothall = projection of search dir along downhill gradient
    // if search direction is not downhill, exit with error

    fdothme = 0.0;
    for (i = 0; i < nvec; i++) fdothme += fvec[i] * h[i];
    if (fdothme <= 0.0) return DOWNHILL;

    // set alphamax so no dof is changed by more than max allowed amount
    // for atom coords, max amount = dmax
    // for extra per-atom dof, max amount = extra_max[]
    // for extra global dof, max amount is set by fix
    // also ensure alphamax <= ALPHA_MAX
    // else will have to backtrack from huge value when forces are tiny
    // if all search dir components are already 0.0, exit with error

    hme = 0.0;
    for (i = 0; i < nvec; i++) hme = max(hme, fabs(h[i]));
    alphamax = std::min(ALPHA_MAX, dist_max/hme);
    if (hme == 0.0) return ZEROFORCE;

    // store box and values of all dof at start of linesearch

    for (i = 0; i < nvec; i++) x0[i] = xvec[i];

    // backtrack with alpha until energy decrease is sufficient
    // or until get to small energy change, then perform quadratic projection

    alpha = alphamax;
    fhprev = fdothme;
    engprev = orig_ene;
    alphaprev = 0.0;

    // // important diagnostic: test the gradient against energy
    // prec_t etmp;
    // prec_t alphatmp = alphamax*1.0e-4;
    // etmp = alpha_step(alphatmp,1);
    // printf("alpha = %g dele = %g dele_force = %g err = %g\n",
    //        alphatmp,etmp-eoriginal,-alphatmp*fdothall,
    //        etmp-eoriginal+alphatmp*fdothall);
    // alpha_step(0.0,1);

    while (true) {
        curr_ene = alpha_step(alpha);

        // compute new fh, alpha, delfh

        dot[0] = dot[1] = 0.0;
        for (i = 0; i < nvec; i++)
        {
            dot[0] += fvec[i] * fvec[i] * mvec[i];
            dot[1] += fvec[i] * h[i];
        }
        fh = dot[1];
        delfh = fh - fhprev;

        // if fh or delfh is epsilon, reset to starting point, exit with error

        if (fabs(fh) < EPS_QUAD || fabs(delfh) < EPS_QUAD)
        {
            curr_ene = alpha_step(0.0);
            return ZEROQUAD;
        }

        // Check if ready for quadratic projection, equivalent to secant method
        // alpha0 = projected alpha

        relerr = fabs(1.0-(0.5*(alpha-alphaprev)*(fh+fhprev)+curr_ene)/engprev);
        alpha0 = alpha - (alpha-alphaprev)*fh/delfh;

        if (relerr <= QUADRATIC_TOL && alpha0 > 0.0 && alpha0 < alphamax)
        {
            curr_ene = alpha_step(alpha0);
            if (curr_ene - orig_ene < EMACH) {
                return 0;
            }
        }

        // if backtracking energy change is better than ideal, exit with success

        de_ideal = -BACKTRACK_SLOPE*alpha*fdothme;
        de = curr_ene - orig_ene;

        if (de <= de_ideal) return 0;

        // save previous state

        fhprev = fh;
        engprev = curr_ene;
        alphaprev = alpha;

        // reduce alpha

        alpha *= ALPHA_REDUCE;

        // backtracked all the way to 0.0
        // reset to starting point, exit with error

        if (alpha <= 0.0 || de_ideal >= -EMACH) {
            curr_ene = alpha_step(0.0);
            return ZEROALPHA;
        }
    }
}

/* ----------------------------------------------------------------------

forcezero linesearch method - seeks a zero of force in a robust manner.
    (motivated by a line minimization routine of f77 DYNAMO code)

central idea:
  In each linesearch we attempt to converge to a zero of force
  (usual case) or reduces forces (worst case).
  Energy does not play any role in the search procedure,
  except we ensure that it doesn't increase.

pseudo code:
  i)  Fix an alpha max:
        // also account for nextra atom & global
        alpha_max <= dmax/hmaxall

  ii) Initialize:
        fhCurr = current_force.dot.search_direction
        fhoriginal = fhCurr
        // try decreasing the energy to 1/10 of initial
        alpha_init = 0.1*fabs(eoriginal)/fhCurr;
        // initial alpha is smaller than alpha_max
        alpha_del = MIN(alpha_init, 0.5*alpha_max);
        alpha = 0.0
  iii) Loop:
        backtrack = false
        alpha += alpha_del
        if (alpha > alpha_max):
           // we have done enough in the search space
           EXIT with success

        Step with the new alpha
        Compute:
           current energy and 'fhCurr'
           de = ecurrent - eprev

           // ZERO_ENERGY = 1e-12, is max allowed energy increase
           if (de > ZERO_ENERGY):
              bactrack = true

           // GRAD_TOL = 0.1
           if ((not backtrack) && (fabs(fhCurr/fh0) <= GRAD_TOL)):
              // forces sufficiently reduced without energy increase
              EXIT with success

           // projected force changed sign but didn't become small enough
           if ( fhCurr < 0):
              backtrack = true

           if (bactrack):
              // forces along search direction changed sign
              if (fhCurr < 0):
                 Get alpha_del by solving for zero
                    of force (1D Newton's Method)
              else:
                 // force didn't change sign but only energy increased,
                 // we overshot a minimum which is very close to a
                 // maximum (or there is an inflection point)

                 // New alpha_del should be much smaller
                 // ALPHA_FACT = 0.1
                 alpha_del *= ALPHA_FACT

                 // Check to see if new 'alpha_del' isn't too small
                 if (alpha_del < MIN_ALPHA):
                    EXIT with failure("linesearch alpha is zero")

               Undo the step of alpha.

           // continue the loop with a new alpha_del
           else:
              Get new alpha_del by linearizing force and solving for its zero

 ---------------------------------------------------------------------- */

index_t KMC::minLineSearch::linemin_forcezero(prec_t orig_ene, prec_t &alpha)
{
    int i,m,n;
    prec_t fdothme,hme;
    prec_t de;

    prec_t alpha_max, alpha_init, alpha_del;
    // projection of: force on itself, current force on search direction,
    prec_t ffCurr, fhCurr;
    // previous force on search direction, initial force on search direction
    prec_t fhPrev, fhoriginal;
    // current energy, previous energy
    prec_t engCurr, engPrev;
    bool backtrack;

    // hardcoded constants

    // factor by which alpha is reduced when backtracking
    prec_t ALPHA_FACT = 0.1;
    // maximum amount by which we'll permit energy increase
    prec_t ZERO_ENERGY = 1e-12;
    // fraction to which we want to reduce the directional derivative
    prec_t GRAD_TOL = 0.1;
    // largest alpha increment which will trigger a failed_linesearch
    prec_t MIN_ALPHA_FAC = 1e-14;
    prec_t LIMIT_BOOST = 4.0;

    // fdothall = projection of search dir along downhill gradient
    // if search direction is not downhill, exit with error

    fdothme = 0.0;
    for (i = 0; i < nvec; i++) fdothme += fvec[i]*h[i];
    if (fdothme <= 0.0) return DOWNHILL;

    // set alpha so no dof is changed by more than max allowed amount
    // for atom coords, max amount = dmax
    // for extra per-atom dof, max amount = extra_max[]
    // for extra global dof, max amount is set by fix

    // also ensure alpha <= ALPHA_MAX else will have
    // to backtrack from huge value when forces are tiny

    // if all search dir components are already 0.0, exit with error

    hme = 0.0;
    for (i = 0; i < nvec; i++) hme = max(hme,fabs(h[i]));
    alpha_max = dist_max/hme;
    if (hme == 0.0) return ZEROFORCE;

    // initialize important variables before main linesearch loop

    ffCurr = 0.0;
    fhCurr = fdothme;
    fhoriginal = fhCurr;
    engCurr = orig_ene;

    // stores energy difference due to the current move

    de = 0.0;

    // choosing the initial alpha that we'll use
    // rough estimate that'll decrease energy to 1/10

    alpha_init = 0.1*fabs(orig_ene)/fdothme;

    // initialize aplha to 0.0

    alpha = 0.0;

    // compute increment to alpha, ensure that we
    // don't take the largest allowed alpha
    // first alpha that will actually apply

    alpha_del = std::min(alpha_init,0.5*alpha_max);

    // main linesearch loop

    while (true) {
        backtrack = false;
        fhPrev = fhCurr;
        engPrev = engCurr;

        // apply the increment to alpha, but first
        // check whether we are still in allowed search space

        alpha += alpha_del;
        if (alpha > alpha_max) {
            // undo the increment

            alpha -= alpha_del;

            // exit linesearch with success: have done
            // enough in allowed search space

            return 0;
        }

        // move the system

        // '1' updates coordinates of atoms which cross PBC

        engCurr = alpha_step(alpha);
        curr_ene = engCurr;

        // compute the new directional derivative and also f_dot_f

        fhCurr = compute_dir_deriv(ffCurr);

        // energy change

        de = engCurr - engPrev;

        // if the function value increases measurably,
        // then we have to reduce alpha

        if (de >= ZERO_ENERGY)
        backtrack = true;

        // check if the directional derivative has sufficiently decreased
        // NOTE: the fabs is essential here

        if ((!backtrack) && (fabs(fhCurr/fhoriginal) <= GRAD_TOL)) return 0;

        // check if the directional derivative changed sign
        // but it's not small: we overshot the minima -- BACKTRACK

        if (fhCurr < 0.0) backtrack = true;

        // backtrack by undoing step and choosing a new alpha

        if (backtrack) {

            // move back

            alpha -= alpha_del;

            // choose new alpha
            // if the force changed sign, linearize force and
            // solve for new alpha_del

            if (fhCurr < 0.0)
                alpha_del *= fhPrev/(fhPrev - fhCurr);
            else

            // force didn't change sign but only energy increased,
            // we overshot a minimum which is very close to a maxima
            // (or there is an inflection point)
            // new alpha_del should be much smaller

                alpha_del *= ALPHA_FACT;

            // since we moved back ...

            engCurr = engPrev;
            curr_ene = engCurr;
            fhCurr = fhPrev;

            // if new move is too small then we have failed;
            // exit with 'failed_linesearch'

            if (hme*alpha_del <= MIN_ALPHA_FAC) {

                // undo all line minization moves

                engCurr = alpha_step(0.0);
                curr_ene = engCurr;
                return ZEROALPHA;
            }

        } else {

            // get a new alpha by linearizing force and start over

            prec_t boostFactor = LIMIT_BOOST;

            // avoids problems near an energy inflection point

            if (fhPrev > fhCurr)
                boostFactor = fhCurr/(fhPrev - fhCurr);

            // don't want to boost too much

            boostFactor = std::min(boostFactor, LIMIT_BOOST);
            alpha_del *= boostFactor;
        }
    }
}

/* ---------------------------------------------------------------------- */

prec_t KMC::minLineSearch::alpha_step(prec_t alpha)
{
    int i;

    // reset to starting point

    for (i = 0; i < nvec; i++) xvec[i] = x0[i];

    // step forward along h

    if (alpha > 0.0)
    {
        for (i = 0; i < nvec; i++) xvec[i] += alpha*h[i];
    }

    // compute and return new energy

    ++neval;
    return energy_force();
}

/* ---------------------------------------------------------------------- */

// compute projection of force on: itself and the search direction

prec_t KMC::minLineSearch::compute_dir_deriv(prec_t &ff)
{
   int i,m,n;
   prec_t *hatom, *fatom;
   prec_t dot[2],dotall[2];
   prec_t fh;

   // compute new fh, alpha, delfh

    dot[0] = dot[1] = 0.0;
    for (i = 0; i < nvec; i++) {
        dot[0] += fvec[i]*fvec[i];
        dot[1] += fvec[i]*h[i];
    }

    ff = dot[0];
    fh = dot[1];
    return fh;
}
