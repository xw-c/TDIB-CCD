set_project("subdiv-ccd")
set_languages("cxxlatest")
set_optimize("faster")

add_rules("mode.debug", "mode.release")

add_requires("eigen","fmt")
add_requires("nlohmann_json")
add_requires("spdlog")

target("rational")
    set_kind("static")
    add_includedirs("rational-cpp/src", {public = true})
    add_includedirs("gmp/include", {public = true})
    add_linkdirs("gmp/lib")
    add_links("libgmp-3")
    add_headerfiles("gmp/include/*.h")
    add_headerfiles("rational-cpp/src/**.hpp")
    add_files      ("rational-cpp/src/**.cpp")

target("ccd_io")
    set_kind("static")
    add_deps("rational")
    add_packages("spdlog", { public = true })
    add_packages("nlohmann_json", { public = true })
    add_includedirs("CCD-Query-IO/src", {public = true})
    add_headerfiles("CCD-Query-IO/src/**.hpp")
    add_files      ("CCD-Query-IO/src/**.cpp")

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

-- target("mesh")
--     set_kind("binary")
--     add_includedirs("core", {public = true})
--     add_packages("eigen", {public = true})
--     -- add_deps("viewer")
--     add_files("SimpleMesh/scene.cpp")

-- target("interval")
--     set_kind("binary")
--     add_includedirs("core", {public = true})
--     add_packages("eigen", {public = true})
--     -- add_deps("viewer")
--     add_files("scene/interval.cpp")

target("benchmark")
    set_kind("binary")
    add_deps("ccd_io")
    add_includedirs("core", {public = true})
    add_packages("eigen","fmt", {public = true})
    add_files("benchmark/**.cpp")