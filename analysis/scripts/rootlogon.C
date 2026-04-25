// rootlogon.C — set up ACLiC include/link paths for prad2 analysis scripts
//
// Two modes, auto-detected:
//
//   1. BUILD-tree mode  (preferred when CMakeCache.txt is reachable)
//        cd <build_dir>
//        root -l ../analysis/scripts/rootlogon.C
//      Or set PRAD2_BUILD_DIR to run from anywhere.  Reads CMakeCache.txt
//      to locate the source tree, the json fetched dep, and the resolved
//      CODA paths; links against the .a archives sitting in the build dir.
//
//   2. INSTALL-tree mode  (when no CMakeCache.txt is found)
//        source <prefix>/bin/prad2_setup.sh        # sets PRAD2_DATABASE_DIR
//        root -l <prefix>/share/prad2evviewer/analysis/scripts/rootlogon.C
//      Derives the install prefix from PRAD2_DATABASE_DIR (the parent of
//      share/prad2evviewer/database), then picks up libs from
//      <prefix>/lib{,64}/ and headers from <prefix>/include/{prad2dec,
//      prad2det,prad2analysis,nlohmann}.  Falls back to the Hall-B CODA
//      path for libevio.a if not installed alongside (override with
//      PRAD2_CODA_ROOT).
//
// Each path probe is logged as "[+] tag : path" (found) or "[-] tag : path"
// (skipped/not present) so a failed setup is easy to debug.  Set
// PRAD2_ROOTLOGON_QUIET=1 to suppress the per-probe lines.

{
    // -------------------------------------------------------------------------
    // Verbose probe helpers — VLOG is a macro (not a generic lambda) so this
    // file works under every cling version since ROOT 6.0.  Set
    // PRAD2_ROOTLOGON_QUIET=1 to suppress the per-probe lines.
    // -------------------------------------------------------------------------
    bool gQuiet = !TString(gSystem->Getenv("PRAD2_ROOTLOGON_QUIET")).IsNull();

    #define VLOG(...)  do { if (!gQuiet) Printf(__VA_ARGS__); } while (0)

    // Probe a single path and print "[+|-] tag : path"; return true if found.
    auto probe = [&](const char *tag, const TString &path) -> bool {
        bool ok = !path.IsNull() && !gSystem->AccessPathName(path);
        VLOG("  [%s] %-14s %s", ok ? "+" : "-", tag, path.Data());
        return ok;
    };
    // Walk the candidate list, log each with [+|-], return the first hit.
    auto pickFirst = [&](const char *tag,
                         std::initializer_list<TString> candidates) -> TString {
        TString winner;
        for (const auto &p : candidates) {
            bool ok = !p.IsNull() && !gSystem->AccessPathName(p);
            VLOG("  [%s] %-14s %s", ok ? "+" : "-", tag, p.Data());
            if (ok && winner.IsNull()) winner = p;
        }
        return winner;
    };

    // -------------------------------------------------------------------------
    // Mode probe
    // -------------------------------------------------------------------------
    TString buildDir = gSystem->Getenv("PRAD2_BUILD_DIR");
    if (buildDir.IsNull()) buildDir = gSystem->pwd();

    Printf("==== prad2 rootlogon ====");
    VLOG("[probe] mode detection");
    bool inBuildMode = probe("CMakeCache",
                             Form("%s/CMakeCache.txt", buildDir.Data()));

    TString sourceDir;        // build mode
    TString prefix;           // install mode
    TString incDec, incDet, incAna, incJson, incCoda;
    TString libDec, libDet, libEvio;
    TString modeLabel = inBuildMode ? "build" : "install";

    if (inBuildMode) {
        // ---------------------------------------------------------------------
        // BUILD MODE
        // ---------------------------------------------------------------------
        Printf("[mode]  build-tree at %s", buildDir.Data());

        // --- source dir from CMakeCache.txt ---
        VLOG("[probe] source dir from CMakeCache.txt");
        std::ifstream cache(Form("%s/CMakeCache.txt", buildDir.Data()));
        std::string line;
        while (std::getline(cache, line)) {
            if (line.find("CMAKE_HOME_DIRECTORY") != std::string::npos ||
                line.find("prad2evviewer_SOURCE_DIR") != std::string::npos) {
                auto eq = line.find('=');
                if (eq != std::string::npos) {
                    sourceDir = line.substr(eq + 1).c_str();
                    break;
                }
            }
        }
        if (sourceDir.IsNull()) {
            VLOG("  [-] CMAKE_HOME_DIRECTORY not in CMakeCache.txt — guessing");
            sourceDir = gSystem->DirName(gSystem->DirName(
                gSystem->Which(".", "analysis/scripts/rootlogon.C")));
            if (sourceDir.IsNull()) sourceDir = "..";
        }
        VLOG("  [+] %-14s %s", "sourceDir", sourceDir.Data());

        // --- CODA paths from CMakeCache.txt ---
        VLOG("[probe] CODA / evio entries in CMakeCache.txt");
        TString codaLibDir, evioLibPath;
        {
            std::ifstream cache2(Form("%s/CMakeCache.txt", buildDir.Data()));
            std::string ln;
            while (std::getline(cache2, ln)) {
                auto eq = ln.find('=');
                if (eq == std::string::npos) continue;
                TString val = ln.substr(eq + 1).c_str();
                if (ln.find("CODA_INCLUDE_DIR:") != std::string::npos) incCoda = val;
                else if (ln.find("CODA_LIB_DIR:")  != std::string::npos) codaLibDir = val;
                else if (ln.find("EVIO_LIB:")       != std::string::npos) evioLibPath = val;
            }
        }
        VLOG("  [%s] %-14s %s",
             incCoda.IsNull()  ? "-" : "+", "CODA_INCLUDE",
             incCoda.IsNull()  ? "(not set)"  : incCoda.Data());
        VLOG("  [%s] %-14s %s",
             codaLibDir.IsNull() ? "-" : "+", "CODA_LIB",
             codaLibDir.IsNull() ? "(not set)"  : codaLibDir.Data());
        VLOG("  [%s] %-14s %s",
             evioLibPath.IsNull() ? "-" : "+", "EVIO_LIB",
             evioLibPath.IsNull() ? "(not set)" : evioLibPath.Data());

        // --- header dirs (build tree) ---
        incDec  = Form("%s/prad2dec/include",  sourceDir.Data());
        incDet  = Form("%s/prad2det/include",  sourceDir.Data());
        incAna  = Form("%s/analysis/include",  sourceDir.Data());

        // --- nlohmann/json (FetchContent) ---
        VLOG("[probe] nlohmann/json (FetchContent)");
        incJson = pickFirst("json-include", {
            Form("%s/_deps/json-src/include",        buildDir.Data()),
            Form("%s/_deps/json-src/single_include", buildDir.Data()),
        });

        // --- find archives in the build tree (or fall back to CODA) ---
        VLOG("[probe] libprad2dec.a");
        libDec = pickFirst("prad2dec.a", {
            Form("%s/lib/lib%s.a",       buildDir.Data(), "prad2dec"),
            Form("%s/lib%s.a",           buildDir.Data(), "prad2dec"),
            Form("%s/prad2dec/lib%s.a",  buildDir.Data(), "prad2dec"),
            Form("%s/lib/lib%s.so",      buildDir.Data(), "prad2dec"),
        });
        VLOG("[probe] libprad2det.a");
        libDet = pickFirst("prad2det.a", {
            Form("%s/lib/lib%s.a",       buildDir.Data(), "prad2det"),
            Form("%s/lib%s.a",           buildDir.Data(), "prad2det"),
            Form("%s/prad2det/lib%s.a",  buildDir.Data(), "prad2det"),
            Form("%s/lib/lib%s.so",      buildDir.Data(), "prad2det"),
        });
        VLOG("[probe] libevio.a");
        if (!evioLibPath.IsNull() && probe("evio (cache)", evioLibPath)) {
            libEvio = evioLibPath;
        } else {
            // No cache hint or it didn't resolve — sweep build dir, then
            // fall back to the CODA lib dir parsed from CMakeCache.txt.
            std::vector<TString> evCandidates = {
                Form("%s/lib/libevio.a", buildDir.Data()),
                Form("%s/libevio.a",     buildDir.Data()),
            };
            if (!codaLibDir.IsNull())
                evCandidates.push_back(Form("%s/libevio.a", codaLibDir.Data()));
            for (const auto &p : evCandidates) {
                bool ok = !p.IsNull() && !gSystem->AccessPathName(p);
                VLOG("  [%s] %-14s %s", ok ? "+" : "-", "evio.a", p.Data());
                if (ok && libEvio.IsNull()) libEvio = p;
            }
        }
    }
    else {
        // ---------------------------------------------------------------------
        // INSTALL MODE
        // ---------------------------------------------------------------------
        TString dbDir = gSystem->Getenv("PRAD2_DATABASE_DIR");
        if (dbDir.IsNull()) {
            Printf("[ERROR] no CMakeCache.txt at %s and PRAD2_DATABASE_DIR is",
                   buildDir.Data());
            Printf("        not set.  Either cd into your build dir, set");
            Printf("        PRAD2_BUILD_DIR=<build>, or source");
            Printf("        <prefix>/bin/prad2_setup.sh.");
        } else {
            // <prefix>/share/prad2evviewer/database -> <prefix>
            prefix = gSystem->DirName(gSystem->DirName(gSystem->DirName(dbDir)));
            Printf("[mode]  install-tree at %s", prefix.Data());
            VLOG("        derived from PRAD2_DATABASE_DIR=%s", dbDir.Data());

            // Header layout matches install rules in prad2dec/prad2det/analysis
            // CMakeLists.txt; nlohmann/json installs to <prefix>/include/nlohmann
            // automatically via FetchContent.
            incDec  = Form("%s/include/prad2dec",      prefix.Data());
            incDet  = Form("%s/include/prad2det",      prefix.Data());
            incAna  = Form("%s/include/prad2analysis", prefix.Data());
            incJson = Form("%s/include",               prefix.Data());

            // --- libs in <prefix>/lib or lib64 (RHEL convention) ---
            VLOG("[probe] libprad2dec.a");
            libDec = pickFirst("prad2dec.a", {
                Form("%s/lib/libprad2dec.a",   prefix.Data()),
                Form("%s/lib64/libprad2dec.a", prefix.Data()),
            });
            VLOG("[probe] libprad2det.a");
            libDet = pickFirst("prad2det.a", {
                Form("%s/lib/libprad2det.a",   prefix.Data()),
                Form("%s/lib64/libprad2det.a", prefix.Data()),
            });

            // libevio: try install prefix first, then Hall-B CODA default
            // (override with PRAD2_CODA_ROOT for non-RHEL9 nodes).
            VLOG("[probe] libevio.a");
            libEvio = pickFirst("evio.a", {
                Form("%s/lib/libevio.a",   prefix.Data()),
                Form("%s/lib64/libevio.a", prefix.Data()),
            });
            if (libEvio.IsNull()) {
                const char *codaRoot = gSystem->Getenv("PRAD2_CODA_ROOT");
                if (!codaRoot || !*codaRoot)
                    codaRoot = "/usr/clas12/release/2.0.0/coda/Linux_x86_64_RHEL9";
                VLOG("[probe] libevio.a (CODA fallback at %s)", codaRoot);
                libEvio = pickFirst("evio.a", {
                    Form("%s/lib/libevio.a", codaRoot),
                });
            }
        }
    }

    // -------------------------------------------------------------------------
    // Apply include paths (verbose).
    // -------------------------------------------------------------------------
    Printf("[probe] include paths");
    auto addInc = [&](const char *tag, const TString &p) {
        if (p.IsNull()) {
            VLOG("  [-] %-14s (not resolved)", tag);
            return;
        }
        bool exists = !gSystem->AccessPathName(p);
        VLOG("  [%s] %-14s %s", exists ? "+" : "-", tag, p.Data());
        if (exists) gSystem->AddIncludePath(Form("-I%s", p.Data()));
    };
    addInc("prad2dec",      incDec);
    addInc("prad2det",      incDet);
    addInc("prad2analysis", incAna);
    addInc("nlohmann/json", incJson);
    addInc("CODA",          incCoda);
    if (inBuildMode) {
        TString srcInc = Form("%s/src", sourceDir.Data());
        addInc("server src",  srcInc);
    }

    // -------------------------------------------------------------------------
    // Apply linker line.
    // -------------------------------------------------------------------------
    if (libDec.IsNull() || libDet.IsNull()) {
        Printf("\n[WARN] could not find libprad2dec.a or libprad2det.a.");
        if (inBuildMode)
            Printf("       Build first:  cmake --build %s", buildDir.Data());
        else
            Printf("       Reinstall the project so libraries land in <prefix>/lib.");
        Printf("       Scripts that require ACLiC (.C+) will fail to link.\n");
    } else {
        TString linkLibs = Form("%s %s", libDet.Data(), libDec.Data());
        if (!libEvio.IsNull()) linkLibs += Form(" %s", libEvio.Data());
        linkLibs += " -lexpat";   // evio dependency
        gSystem->AddLinkedLibs(linkLibs);
        Printf("[link]  %s", linkLibs.Data());
    }

    // -------------------------------------------------------------------------
    // Database path — honour an explicit override; otherwise fall back to
    // the build-tree database/ in build mode (install mode already had it).
    // -------------------------------------------------------------------------
    TString dbDir = gSystem->Getenv("PRAD2_DATABASE_DIR");
    if (dbDir.IsNull()) {
        if (inBuildMode) dbDir = Form("%s/database", sourceDir.Data());
        else             dbDir = Form("%s/share/prad2evviewer/database", prefix.Data());
        gSystem->Setenv("PRAD2_DATABASE_DIR", dbDir);
        VLOG("[db]    PRAD2_DATABASE_DIR not set; defaulting to %s", dbDir.Data());
    }
    Printf("[db]    %s", dbDir.Data());

    Printf("==== ready (%s mode) ====", modeLabel.Data());
    Printf("Example: .x gem_clusters_to_root.C+(\"/data/run.evio\", \"out.root\")");
}
