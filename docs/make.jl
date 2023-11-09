using Documenter
using QMOC

DocMeta.setdocmeta!(QMOC, :DocTestSetup, :(using QMOC); recursive=true)

makedocs(
    sitename="QMOC",
    modules=[QMOC]
    authors = "Daniel Simm",
    repo="https://github.com/danielsimm/QMOC.jl/blob/{commit}{path}#{line}",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://danielsimm.github.io/QMOC.jl",
        edit_link="main",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/danielsimm/QMOC.jl",
    devbranch="main"
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
