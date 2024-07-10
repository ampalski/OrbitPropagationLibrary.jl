export OplIn

abstract type OplIn end
abstract type InitialState end

struct Cart_InitialState <: InitialState
    cartState::SVector{6,Float64}
end

struct TwoBody_OplIn <: OplIn
    state0::InitialState
end
