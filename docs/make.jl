using VicsekModel
using Documenter

DocMeta.setdocmeta!(VicsekModel, :DocTestSetup, :(using VicsekModel); recursive=true)

makedocs(;
    modules=[VicsekModel],
    authors="Masashi Shiraishi <shiraishi.mu@gmail.com> and contributors",
    sitename="VicsekModel.jl",
    format=Documenter.HTML(;
        canonical="https://shiraishi-mu.github.io/VicsekModel.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/shiraishi-mu/VicsekModel.jl",
    devbranch="main",
)
