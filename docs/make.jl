using Documenter
using QMOC

makedocs(
    sitename = "QMOC",
    format = Documenter.HTML(),
    modules = [QMOC]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
