# test architecture

We use [munit] for C unit testing.

one suite per function, several tests per suite for different mat/vec sizes.
after the first tests for matrix creation, we create setup and teardown functions
to create and destroy mat and vec structs.
suite naming: /blas/struct/<function name>
test naming: /blas/struct/<function name>/<mat or vec size>
