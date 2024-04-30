set_project("subdiv-ccd")
set_languages("cxxlatest")
set_optimize("faster")

add_rules("mode.debug", "mode.release")
add_requires("eigen","fmt")

-- target("demo")
--     set_kind("binary")
--     -- add_includedirs("core", {public = true})
--     add_packages("eigen", {public = true})
--     add_files("backup/ccd/test.cpp")

-- target("core")
--     set_kind("static")
--     add_packages("eigen", {public = true})
--     add_includedirs("core", {public = true})
--     add_files("core/*.cpp")

target("scene")
    set_kind("binary")
    add_includedirs("core", {public = true})
    add_packages("eigen","fmt", {public = true})
    add_files("scene/demo.cpp","core/*.cpp")

target("linear")
    set_kind("binary")
    add_packages("eigen","fmt", {public = true})
    add_files("LinearMesh/demo.cpp", "core/argsParser.cpp")
