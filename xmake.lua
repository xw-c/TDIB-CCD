set_project("subdiv-ccd")
set_languages("cxxlatest")

add_rules("mode.debug", "mode.release")

add_requires("eigen")

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
    add_files("ccd/**.cpp")