// rootlogon.C — set up ACLiC include/link paths for prad2 analysis scripts
//
// Usage:
//   cd <build_dir>                     # must be run from the build directory
//   root -l ../analysis/scripts/rootlogon.C
//   .x ../analysis/scripts/analysis.C+("recon.root")
//
// Or set PRAD2_BUILD_DIR env var to run from anywhere:
//   export PRAD2_BUILD_DIR=/path/to/build
//   root -l analysis/scripts/rootlogon.C
//
// The script auto-detects the source and build directories and configures
// ACLiC to compile .C macros against libprad2dec and libprad2det.

{
    // --- locate build and source directories ---
    TString buildDir = gSystem->Getenv("PRAD2_BUILD_DIR");
    if (buildDir.IsNull()) buildDir = gSystem->pwd();

    // find source dir: look for prad2dec/include relative to CMakeCache
    TString sourceDir;
    std::ifstream cache(Form("%s/CMakeCache.txt", buildDir.Data()));
    if (cache.is_open()) {
        std::string line;
        while (std::getline(cache, line)) {
            if (line.find("CMAKE_HOME_DIRECTORY") != std::string::npos ||
                line.find("prad2evviewer_SOURCE_DIR") != std::string::npos) {
                auto eq = line.find('=');
                if (eq != std::string::npos)
                    sourceDir = line.substr(eq + 1).c_str();
                break;
            }
        }
    }
    if (sourceDir.IsNull()) {
        // fallback: assume build is a subdirectory of source
        sourceDir = gSystem->DirName(gSystem->DirName(
            gSystem->Which(".", "analysis/scripts/rootlogon.C")));
        if (sourceDir.IsNull()) sourceDir = "..";
    }

    Printf("PRad2 source : %s", sourceDir.Data());
    Printf("PRad2 build  : %s", buildDir.Data());

    // --- parse CMakeCache.txt for CODA and json paths -------------------------
    TString codaIncludeDir, codaLibDir, evioLibPath;
    {
        std::ifstream cache2(Form("%s/CMakeCache.txt", buildDir.Data()));
        std::string ln;
        while (std::getline(cache2, ln)) {
            auto eq = ln.find('=');
            if (eq == std::string::npos) continue;
            TString val = ln.substr(eq + 1).c_str();
            if (ln.find("CODA_INCLUDE_DIR:") != std::string::npos) codaIncludeDir = val;
            else if (ln.find("CODA_LIB_DIR:")  != std::string::npos) codaLibDir = val;
            else if (ln.find("EVIO_LIB:")       != std::string::npos) evioLibPath = val;
        }
    }

    // --- find nlohmann/json (fetched by cmake) ---
    TString jsonInclude;
    // try common FetchContent paths
    TString jsonPath = Form("%s/_deps/json-src/include", buildDir.Data());
    if (gSystem->AccessPathName(jsonPath)) // returns true if NOT accessible
        jsonPath = Form("%s/_deps/json-src/single_include", buildDir.Data());
    if (!gSystem->AccessPathName(jsonPath))
        jsonInclude = Form("-I%s", jsonPath.Data());

    // --- set include paths ---
    gSystem->AddIncludePath(Form("-I%s/prad2dec/include", sourceDir.Data()));
    gSystem->AddIncludePath(Form("-I%s/prad2det/include", sourceDir.Data()));
    gSystem->AddIncludePath(Form("-I%s/analysis/include", sourceDir.Data()));
    gSystem->AddIncludePath(Form("-I%s/src", sourceDir.Data()));
    if (!jsonInclude.IsNull())
        gSystem->AddIncludePath(jsonInclude);
    if (!codaIncludeDir.IsNull() && !gSystem->AccessPathName(codaIncludeDir))
        gSystem->AddIncludePath(Form("-I%s", codaIncludeDir.Data()));

    // --- find static libraries ---
    // libraries may be in lib/, or directly in subdirectories
    auto findLib = [&](const char *name) -> TString {
        std::vector<TString> candidates = {
            Form("%s/lib/lib%s.a", buildDir.Data(), name),
            Form("%s/lib%s.a", buildDir.Data(), name),
            Form("%s/prad2dec/lib%s.a", buildDir.Data(), name),
            Form("%s/prad2det/lib%s.a", buildDir.Data(), name),
            Form("%s/lib/lib%s.so", buildDir.Data(), name),
        };
        // also check CODA lib dir (Hall-B installation)
        if (!codaLibDir.IsNull())
            candidates.push_back(Form("%s/lib%s.a", codaLibDir.Data(), name));
        for (auto &p : candidates)
            if (!gSystem->AccessPathName(p)) return p;
        return "";
    };

    TString libDec  = findLib("prad2dec");
    TString libDet  = findLib("prad2det");
    // prefer the exact evio path cmake resolved, fall back to search
    TString libEvio;
    if (!evioLibPath.IsNull() && !gSystem->AccessPathName(evioLibPath))
        libEvio = evioLibPath;
    else
        libEvio = findLib("evio");

    if (libDec.IsNull() || libDet.IsNull()) {
        Printf("\n  WARNING: could not find libprad2dec.a or libprad2det.a in %s",
               buildDir.Data());
        Printf("  Make sure you built the project first: cmake --build %s",
               buildDir.Data());
        Printf("  Scripts that require ACLiC (.C+) will fail to link.\n");
    } else {
        TString linkLibs = Form("%s %s", libDet.Data(), libDec.Data());
        if (!libEvio.IsNull()) linkLibs += Form(" %s", libEvio.Data());
        linkLibs += " -lexpat";  // evio dependency

        gSystem->AddLinkedLibs(linkLibs);
        Printf("ACLiC libs   : %s", linkLibs.Data());
    }

    // --- set database path ---
    TString dbDir = gSystem->Getenv("PRAD2_DATABASE_DIR");
    if (dbDir.IsNull()) dbDir = Form("%s/database", sourceDir.Data());
    gSystem->Setenv("PRAD2_DATABASE_DIR", dbDir);
    Printf("Database     : %s", dbDir.Data());

    Printf("\nReady. Examples:");
    Printf("  .x %s/analysis/scripts/lms_alpha_normalize.C+(\"/data/run\", 1234)", sourceDir.Data());
}
