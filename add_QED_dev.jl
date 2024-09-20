using Pkg

Pkg.add(; url="https://github.com/QEDjl-project/QEDbase.jl", rev="dev")
@warn "This repository depends on the dev branch of QEDbase.jl\n It is NOT ready for release!"

Pkg.add(; url="https://github.com/QEDjl-project/QEDprocesses.jl.git", rev="dev")
@warn "This repository depends on the dev branch of QEDprocesses.jl\n It is NOT ready for release!"
