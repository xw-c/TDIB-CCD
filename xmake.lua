set_project("subdiv-ccd")
set_languages("cxxlatest")

add_rules("mode.debug", "mode.release")

add_requires("eigen")

target("core")
    set_kind("static")
    -- add_files("core/*.cpp")
    add_includedirs("core", {public = true})
    add_packages("eigen", {public = true})

target("test")
    set_kind("binary")
    add_packages("eigen")

    add_includedirs("dcd")
    add_headerfiles("dcd/**.h")
    add_files("dcd/**.cpp")

target("ccd")
    set_kind("binary")
    add_packages("eigen")

    add_includedirs("ccd")
    add_headerfiles("ccd/**.h")
    add_files("ccd/test.cpp")
    add_deps("core")

target("demo")
    set_kind("binary")
    add_files("ccd/testtri.cpp")
    add_deps("core")