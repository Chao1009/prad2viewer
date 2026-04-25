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
//      prad2det,prad2ana,nlohmann}.  Falls back to the Hall-B CODA
//      path for libevio.a if not installed alongside (override with
//      PRAD2_CODA_ROOT).
//
// Each path probe is logged as "[+] tag : path" (found) or "[-] tag : path"
// (skipped/not present) so a failed setup is easy to debug.
//
// Environment overrides (in order of priority):
//   PRAD2_BUILD_DIR        — build dir if not the cwd
//   PRAD2_DATABASE_DIR     — install-mode prefix anchor (set by prad2_setup.sh)
//   PRAD2_EVIO_LIB         — explicit libevio.a path (skips all evio probes)
//   PRAD2_CODA_ROOT        — non-default CODA install root
//   PRAD2_ROOTLOGON_QUIET  — suppress the per-probe lines

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
    TString libDec, libDet, libAna, libEvio;
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
        VLOG("[probe] libprad2ana.a");
        libAna = pickFirst("prad2ana.a", {
            Form("%s/lib/lib%s.a",       buildDir.Data(), "prad2ana"),
            Form("%s/lib%s.a",           buildDir.Data(), "prad2ana"),
            Form("%s/analysis/lib%s.a",  buildDir.Data(), "prad2ana"),
            Form("%s/lib/lib%s.so",      buildDir.Data(), "prad2ana"),
        });
        // libevio.a — same Hall-B-first-then-fetch shape that
        // prad2dec/CMakeLists.txt uses at configure time, repeated here so
        // ACLiC can resolve the link target whichever path the build took:
        //
        //   1. PRAD2_EVIO_LIB env var (explicit override)
        //   2. EVIO_LIB from CMakeCache.txt (cmake's own resolved path)
        //   3. CMakeCache CODA_LIB_DIR (Hall-B install)
        //   4. fetched evio under <build>/_deps/evio-build/ (FetchContent)
        //   5. <build>/lib/ or <build>/ (other install layouts)
        //
        // The fetched evio's exact archive path varies between releases, so
        // we enumerate the layouts we've seen.  If your build lands it
        // somewhere else, point PRAD2_EVIO_LIB at it and skip the search.
        VLOG("[probe] libevio.a");
        const char *envEvio = gSystem->Getenv("PRAD2_EVIO_LIB");
        if (envEvio && *envEvio && probe("evio (env)", envEvio)) {
            libEvio = envEvio;
        }
        else if (!evioLibPath.IsNull() && probe("evio (cache)", evioLibPath)) {
            libEvio = evioLibPath;
        } else {
            std::vector<TString> evCandidates;
            if (!codaLibDir.IsNull())
                evCandidates.push_back(Form("%s/libevio.a", codaLibDir.Data()));
            // FetchContent layouts (evio-6.x).  The legacy Hall-B-style
            // build puts archives under src/libsrc/<arch>/; the modernized
            // CMake build puts them at the top of evio-build/lib/ or just
            // evio-build/.
            evCandidates.push_back(Form("%s/_deps/evio-build/lib/libevio.a",
                                        buildDir.Data()));
            evCandidates.push_back(Form("%s/_deps/evio-build/libevio.a",
                                        buildDir.Data()));
            evCandidates.push_back(Form("%s/_deps/evio-build/src/libsrc/libevio.a",
                                        buildDir.Data()));
            evCandidates.push_back(Form("%s/_deps/evio-build/src/libsrc/Linux-x86_64/libevio.a",
                                        buildDir.Data()));
            evCandidates.push_back(Form("%s/lib/libevio.a", buildDir.Data()));
            evCandidates.push_back(Form("%s/libevio.a",     buildDir.Data()));
            for (const auto &p : evCandidates) {
                bool ok = !p.IsNull() && !gSystem->AccessPathName(p);
                VLOG("  [%s] %-14s %s", ok ? "+" : "-", "evio.a", p.Data());
                if (ok && libEvio.IsNull()) libEvio = p;
            }
            // Last resort — recursive scan under _deps/evio-build/ for any
            // libevio.a we missed.  Bounded depth, executed only when the
            // enumerated candidates all failed, so this is cheap.
            if (libEvio.IsNull()) {
                TString depsDir = Form("%s/_deps/evio-build", buildDir.Data());
                if (!gSystem->AccessPathName(depsDir)) {
                    VLOG("  [.] %-14s scanning %s ...", "evio.a", depsDir.Data());
                    TString cmd = Form("find '%s' -maxdepth 5 -name libevio.a "
                                       "-print -quit 2>/dev/null", depsDir.Data());
                    TString hit = gSystem->GetFromPipe(cmd);
                    hit = hit.Strip(TString::kBoth);
                    if (!hit.IsNull() && !gSystem->AccessPathName(hit)) {
                        VLOG("  [+] %-14s %s", "evio.a", hit.Data());
                        libEvio = hit;
                    } else {
                        VLOG("  [-] %-14s (no libevio.a under _deps/evio-build)",
                             "evio.a");
                    }
                }
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
            incAna  = Form("%s/include/prad2ana",      prefix.Data());
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
            VLOG("[probe] libprad2ana.a");
            libAna = pickFirst("prad2ana.a", {
                Form("%s/lib/libprad2ana.a",   prefix.Data()),
                Form("%s/lib64/libprad2ana.a", prefix.Data()),
            });

            // libevio: PRAD2_EVIO_LIB override > install prefix (covers both
            // FetchContent-bundled evio and CODA-installed) > Hall-B CODA
            // default at /usr/clas12/.../RHEL9 (override with PRAD2_CODA_ROOT
            // for other arch/distro combos).
            VLOG("[probe] libevio.a");
            const char *envEvio = gSystem->Getenv("PRAD2_EVIO_LIB");
            if (envEvio && *envEvio && probe("evio (env)", envEvio)) {
                libEvio = envEvio;
            } else {
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
    addInc("prad2ana",      incAna);
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
        // Order matters for static-archive linking: a symbol referenced by
        // libprad2ana lives in libprad2det/dec, so analysis comes first
        // (left), det/dec to the right.  evio + expat are leaf deps.
        TString linkLibs;
        if (!libAna.IsNull()) linkLibs += Form("%s ", libAna.Data());
        linkLibs += Form("%s %s", libDet.Data(), libDec.Data());
        if (!libEvio.IsNull()) linkLibs += Form(" %s", libEvio.Data());
        linkLibs += " -lexpat";   // evio dependency
        gSystem->AddLinkedLibs(linkLibs);
        Printf("[link]  %s", linkLibs.Data());
        if (libAna.IsNull())
            Printf("[note]  libprad2ana.a not found — scripts that call "
                   "analysis::* (PhysicsTools / MatchingTools / Replay) will "
                   "fail to link.  Build the analysis target.");
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
