#pragma once
#ifndef FIXRLCUT_H
#define FIXRLCUT_H

#include "FixPair.h"
#include "xml_func.h"
//! Make FixLJCut available to the pair base class in boost
class EvaluatorWrapper;
void export_FixRLCut();

//! Fix for truncated 3-body interactions
/*!
 * Fix to calculate 3-body interactions of particles. The 3-body potential
 * is defined as
 * \f[
 * V(r_{ij}) = 4 \varepsilon \left[ \left(\frac{\sigma}{r_{ij}}\right)^{12} -
 *                               \left(\frac{\sigma}{r_{ij}}\right)^{6}\right],
 * \f]
 * where \f$ r \f$ is the distance between two particles and \f$ \varepsilon \f$
 * and \f$ \sigma \f$ are the two relevant parameters. The LJ pair interaction
 * is only calculated for particles closer than \$r_{\text{cut}}\$.
 *
 * From the potential, the force can be derived as
 * \f[
 * F(r_{ij}) = 24 \varepsilon \frac{1}{r} \left[ 
 *                          2 \left(\frac{\sigma}{r_{ij}}\right)^{12} -
 *                            \left(\frac{\sigma}{r_{ij}}\right)^{6}
 *                          \right].
 * \f]
 * If \f$F(r_{ij}) < 0\f$, then the force is attractive. Otherwise, it is
 * repulsive.
*/

extern const std::string RLCutType;
class FixRLCut : public FixPair {
    public:
        //! Constructor
        FixRLCut(SHARED(State), std::string handle);

        //! Compute forces
        void compute(int);

        //! Compute single point energy
        void singlePointEng(float *);

        //! Prepare Fix
        /*!
         * \returns Always returns True
         * 
         * This function needs to be called before simulation run.
         */
        bool prepareForRun();

        //! Run after simulation
        /*!
         * This function needs to be called after simulation run.
         */
        bool postRun();

        //! Create restart string
        /*!
         * \param format Format of the pair parameters.
         *
                  * \returns restart chunk string.
         */
        std::string restartChunk(std::string format);


        //! Add new type of atoms
        /*!
         * \param handle Not used
         *
         * This function adds a new particle type to the fix.
         */
        void addSpecies(std::string handle);

        //! Return list of cutoff values
        std::vector<float> getRCuts();
    public:
        void setEvalWrapper();
        void setEvalWrapperOrig();

        const std::string epsHandle; //!< Handle for parameter epsilon
        const std::string sigHandle; //!< Handle for parameter sigma
        const std::string rCutHandle; //!< Handle for parameter rCut
        const std::string alphaHandle; //!< Handle for parameter alpha
        const std::string r_aHandle; //!< Handle for r_a
        const std::string a_gHandle; //!< Handle for a_g
        const std::string sig_gHandle; //!< Handle for sigma_g
        const std::string r_gHandle; //!< Handle for r_g
        std::vector<float> epsilons; //!< vector storing epsilon values
        std::vector<float> sigmas; //!< vector storing sigma values
        std::vector<float> rCuts; //!< vector storing cutoff distance values
        std::vector<float> alphas; //!< vector storing alpha
        std::vector<float> r_as; //!<vector storing r_a
        std::vector<float> a_gs; //!< vector storing a_g
        std::vector<float> sig_gs; //!< vector storing sigma_g
        std::vector<float> r_gs; //!< vector storing r_g

        //EvaluatorLJ evaluator; //!< Evaluator for generic pair interactions

};

#endif
