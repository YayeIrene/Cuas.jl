t =Template(;
    user="YayeIrene",
    authors="Irene Ndindabahizi",
    plugins=[
        License(; name="MIT"),
        Codecov(),
        TravisCI(),
        Coveralls(),
        AppVeyor(),
    ],
)

generate("ErrorBudget", t)
