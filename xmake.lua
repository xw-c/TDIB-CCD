set_project("subdiv-ccd")
set_languages("cxxlatest")

add_rules("mode.debug", "mode.release")

add_requires("eigen")

target("core")
    set_kind("static")
    add_files("core/*.cpp")
    add_includedirs("core", {public = true})
    add_packages("eigen", {public = true})

target("demo")
    set_kind("binary")
    add_files("ccd/testtri.cpp")
    add_deps("core")