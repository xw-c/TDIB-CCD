set_project("subdiv-ccd")
set_languages("cxxlatest")

add_rules("mode.debug", "mode.release")

add_requires("eigen")

target("test")
    set_kind("binary")
    add_packages("eigen")

    add_includedirs("src")
    add_headerfiles("src/**.h")
    add_files("src/**.cpp")
