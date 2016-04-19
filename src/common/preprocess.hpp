#ifndef __PREPROCESS__
#define __PREPROCESS__

#include "common.hpp"
#include "CFID.hpp"
#include "Options.hpp"

namespace tarquin 
{

	class Workspace;
	class CBoswell;
	class CFID;
	class Options;
    void AutoPhase(const coord& proc_coord, CFID& fid, const Options& opts, CBoswell& log, bool water_file=false);

    void AutoPhaseNew(const coord& proc_coord, CFID& fid, const Options& opts, CBoswell& log, bool water_file=false);

    bool AutoReferenceCorr(const coord& proc_coord, Options& options, CFID& fid, bool dyn_mode, bool add_h2o, bool skip_beta_guess, CBoswell& bos_log);

	treal ComputeWaterNormalisation(const coord& proc_coord, const CFID& fid, CBoswell& log);

    treal GetTimeDomainAmplitude(const cvm::cvector& y, const cvm::rvector& t);

/*!
 * Do preprocessing on FID, i.e. water removal, simple phasing, referencing, eddy current, etc.
 */
class Preprocessor
{
	public:

		Preprocessor(Workspace& workspace, CBoswell& log);

		void operator() ();

	private:

		Workspace& m_workspace;

		CBoswell& m_log;
};

} // namespace tarquin

#endif
