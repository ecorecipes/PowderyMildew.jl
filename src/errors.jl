"""
    PowderyMildewError

Abstract supertype of the errors thrown by PowderyMildew.jl. Invalid arguments throw
`ArgumentError` instead.
"""
abstract type PowderyMildewError <: Exception end

"""
    SolverError(msg)

An ODE integration did not finish successfully.
"""
struct SolverError <: PowderyMildewError
    msg::String
end

Base.showerror(io::IO, e::SolverError) = print(io, "SolverError: ", e.msg)
