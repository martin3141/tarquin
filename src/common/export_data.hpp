#ifndef __EXPORTDATA__
#define __EXPORTDATA__

#include "CFID.hpp"
#include "CBasis.hpp"
#include "Workspace.hpp"
#include <string>

using namespace tarquin;

void ExportCsvSpectrum(const std::string& strFilename, const Workspace& workspace);

void ExportCsvSpectraAligned(const std::string& strFilename, const Workspace& workspace, bool mag = false);

void ExportCsvFit(const std::string& strFilename, const Workspace& workspace);

void ExportCsvResults(const std::string& strFilename, const Workspace& workspace);

void ExportTxtResults(const std::string& strFilename, const Workspace& workspace);

void ExportPdfResults(const std::string& strFilename, const Workspace& workspace, int fit_num = 0);

void GetTable(std::ostringstream& table, const Workspace& workspace);

#endif
