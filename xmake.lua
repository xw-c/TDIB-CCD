set_project("subdiv-ccd")
set_languages("cxxlatest")
set_optimize("faster")

add_rules("mode.debug", "mode.release")
add_requires("eigen","fmt")

target("scene")
    set_kind("binary")
    add_includedirs("core", {public = true})
    add_packages("eigen","fmt", {public = true})
    add_files("scene/demo.cpp","core/*.cpp")

-- target("linear")
--     set_kind("binary")
--     add_packages("eigen","fmt", {public = true})
--     add_files("linear_tri/demo.cpp", "core/argsParser.cpp")
