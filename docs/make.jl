using Pkg
using QEDFeynmanDiagrams

project_path = Base.Filesystem.joinpath(Base.Filesystem.dirname(Base.source_path()), "..")
Pkg.develop(; path=project_path)

using Documenter
using DocumenterInterLinks

# setup interlinks
links = InterLinks(
    "ComputableDAGs" => "https://computabledags.github.io/ComputableDAGs.jl/dev/",
    "QEDbase" => "https://qedjl-project.github.io/QEDbase.jl/dev/",
    "QEDcore" => "https://qedjl-project.github.io/QEDcore.jl/dev/",
    "QuantumElectrodynamics" => "https://qedjl-project.github.io/QuantumElectrodynamics.jl/dev/",
)

# setup examples using Literate.jl
using Literate

literate_paths = [
    (
        Base.Filesystem.joinpath(project_path, "docs/src/examples/trident.jl"),
        Base.Filesystem.joinpath(project_path, "docs/src/examples/"),
    ),
    (
        Base.Filesystem.joinpath(project_path, "docs/src/examples/compton.jl"),
        Base.Filesystem.joinpath(project_path, "docs/src/examples/"),
    ),
]

for (file, output_dir) in literate_paths
    Literate.markdown(file, output_dir; documenter=true)
    Literate.notebook(file, output_dir)
end

pages = [
    "index.md",
    "Manual" => "manual.md",
    "Examples" => ["Trident" => "examples/trident.md", "Compton" => "examples/compton.md"],
    "Library" => ["Public" => "lib/public.md", "Internal" => "lib/internal.md"],
    "Contribution" => "contribution.md",
]

makedocs(;
    modules=[QEDFeynmanDiagrams],
    checkdocs=:exports,
    authors="Anton Reinhard",
    repo=Documenter.Remotes.GitHub("QEDjl-project", "QEDFeynmanDiagrams.jl"),
    sitename="QEDFeynmanDiagrams.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://QEDjl-project.github.io/QEDFeynmanDiagrams.jl",
        assets=String[],
    ),
    pages=pages,
    plugins=[links],
)
deploydocs(; repo="github.com/QEDjl-project/QEDFeynmanDiagrams.jl.git", push_preview=false)
