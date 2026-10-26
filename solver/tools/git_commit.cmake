# Writes OUT (git_commit.h): the commit of the repository REPO, with "-dirty" when it has uncommitted changes.
# Run at every build; the header is only rewritten when the commit changes, so only what includes it is rebuilt.
find_package(Git QUIET)
set(commit "unknown")
if(GIT_FOUND)
  execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${REPO}" rev-parse --short=12 HEAD
                  OUTPUT_VARIABLE head RESULT_VARIABLE failed OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  if(NOT failed)
    set(commit "${head}")
    execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${REPO}" status --porcelain --untracked-files=no
                    OUTPUT_VARIABLE changes OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if(changes)
      set(commit "${commit}-dirty")
    endif()
  endif()
endif()
set(content "#pragma once\n#define pyTEC3D_GIT_COMMIT \"${commit}\"\n")
if(EXISTS "${OUT}")
  file(READ "${OUT}" old)
endif()
if(NOT "${old}" STREQUAL "${content}")
  file(WRITE "${OUT}" "${content}")
endif()
