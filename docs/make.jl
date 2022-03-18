using Documenter, Cuas
push!(LOAD_PATH,"../src/")
makedocs(sitename="My Documentation")

makedocs(
         sitename = "Cuas.jl",
         modules  = [Cuas],
         pages=[
                "Home" => "index.md"
               ])
deploydocs(;
    repo="github.com/YayeIrene/Cuas.jl",
)


