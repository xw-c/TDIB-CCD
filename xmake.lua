set_project("subdiv-ccd")
set_languages("cxxlatest")
set_optimize("faster")

add_rules("mode.debug", "mode.release")

add_requires("eigen")

target("demo")
    set_kind("binary")
    add_includedirs("core", {public = true})
    add_packages("eigen", {public = true})
    add_files("demo.cpp")