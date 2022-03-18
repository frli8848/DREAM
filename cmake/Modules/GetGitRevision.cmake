# Grab the version number from the git repo.
#
# http://brianmilco.blogspot.no/2012/11/cmake-automatically-use-git-tags-as.html
include(GetGitRevisionDescription)

function(GIT_VERSION prefix)
  git_describe(VERSION tag)

  message(STATUS "Generating version info from: " ${VERSION})
  # Parse the version information into pieces (version tags must have the format v0.3.4)
  string(REGEX REPLACE "^v([0-9]+)\\..*" "\\1" major "${VERSION}")
  string(REGEX REPLACE "^v[0-9]+\\.([0-9]+).*" "\\1" minor "${VERSION}")
  string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.((RC)?[0-9]+).*" "\\1" patch "${VERSION}")
  string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.(RC)?[0-9]+-([0-9]+)-.*" "\\2" commit "${VERSION}")
  string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.(RC)?[0-9]+-[0-9]+-(.*)" "\\2" sha1 "${VERSION}")

  set(${prefix}_VERSION_SHORT "${major}.${minor}.${patch}" PARENT_SCOPE)
  set(${prefix}_VERSION_LONG "${major}.${minor}.${patch}.${commit}" PARENT_SCOPE)
  set(${prefix}_VERSION_MAJOR ${major} PARENT_SCOPE)
  set(${prefix}_VERSION_MINOR ${minor} PARENT_SCOPE)
  set(${prefix}_VERSION_COMMIT ${commit} PARENT_SCOPE)
  set(${prefix}_VERSION_PATCH ${patch} PARENT_SCOPE)
  set(${prefix}_VERSION_SHA1 ${sha1} PARENT_SCOPE)
  set(${prefix}_VERSION "${major}.${minor}.${patch}-${commit}-${sha1}" PARENT_SCOPE)
endfunction()
