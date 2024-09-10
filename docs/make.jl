using Pkg

project_path = Base.Filesystem.joinpath(Base.Filesystem.dirname(Base.source_path()), "..")
Pkg.develop(; path=project_path)

using Documenter
using QEDFeynmanDiagrams

pages = [
    "index.md",
    "Manual" => "manual.md",
    "Library" => ["Public" => "lib/public.md", "Internal" => "lib/internal.md"],
    "Contribution" => "contribution.md",
]

makedocs(;
    modules=[QEDFeynmanDiagrams],
    checkdocs=:exports,
    authors="Anton Reinhard",
    repo=Documenter.Remotes.GitHub("QEDFeynmanDiagrams", "QEDFeynmanDiagrams.jl"),
    sitename="QEDFeynmanDiagrams.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://QEDjl-project.github.io/QEDFeynmanDiagrams.jl",
        assets=String[],
    ),
    pages=pages,
)
deploydocs(; repo="github.com/QEDjl-project/QEDFeynmanDiagrams.jl.git", push_preview=false)
