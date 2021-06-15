"""
    function ct_kenya(t,minpoint_y;schoolseffect = 0.,
                    measures_x = 24.,minpoint_x = 55.,returnpoint_x = 254.,schoolsreopen = 319. )

Defines piece-wise linear contact rates.
"""
function ct_kenya(t,minpoint_y;schoolseffect = 0.,
                    measures_x = 24.,minpoint_x = 55.,returnpoint_x = 254.,schoolsreopen = 319. )
    ct = 1.0*(t<measures_x)
    ct += (1.0 - schoolseffect  + ((t-measures_x)*(minpoint_y  - 1.0)/(minpoint_x - measures_x)))*(t ≥ measures_x && t < minpoint_x)
    ct += (minpoint_y - schoolseffect + ((t-minpoint_x)*(1.0   - minpoint_y)/(returnpoint_x-minpoint_x) ))*(t ≥minpoint_x && t < returnpoint_x)
    ct += (1.0 - schoolseffect )*(t ≥ returnpoint_x && t<schoolsreopen)
    ct += 1.0*(t ≥ schoolsreopen)
end

"""
    function twoSEStransmission(du,u,p,t)

In-place creation of vector field for the transmission model
"""
function twoSEStransmission(du,u,p,t)
    #Get state
	S₁,E₁,I₁,R₁,W₁,C₁,S₂,E₂,I₂,R₂,W₂,C₂,C_sero₁,C_sero₂  = u
    #Get parameters
	R₀,α,γ,ϵ,σ,N₁,N₂,ω,schooleffect,ct_min1,ct_min2 = p
    #Adjust for schools being open or shut
    Reff = R₀*(t < 24 || t ≥ 319) + R₀*schooleffect*(t≥ 24 && t< 319) 
    #Rate of visiting POIs for lower and higher SES groups
    ct1 = ct_kenya(t,ct_min1)
    ct2 = ct_kenya(t,ct_min2)
    #Force of infection on lower and higher SES individuals visiting POIs
    λ₁ = γ*Reff*(ct1*ϵ*I₁ + ct2*(1-ϵ)*I₂)
    λ₂ = γ*Reff*(ct1*(1-ϵ)*I₁ + ct2*ϵ*I₂)
    #Dynamics
	du[1] = -ct1*(S₁/N₁)*λ₁
	du[2] = ct1*((S₁+σ*W₁)/N₁)*λ₁ - α*E₁
	du[3] = α*E₁ - γ*I₁
	du[4] = γ*I₁ - ω*R₁
	du[5] = ω*R₁ - σ*ct1*(W₁/N₁)*λ₁
	du[6] = ct1*((S₁+σ*W₁)/N₁)*λ₁
	du[7] = -ct2*(S₂/N₂)*λ₂
	du[8] = ct2*((S₂+σ*W₂)/N₂)*λ₂ - α*E₂
	du[9] = α*E₂ - γ*I₂
	du[10] = γ*I₂ - ω*R₂
	du[11] = ω*R₂ - σ*ct2*(W₂/N₂)*λ₂
	du[12] = ct2*((S₂+σ*W₂)/N₂)*λ₂
    du[13] = ct1*(S₁/N₁)*λ₁
    du[14] = ct2*(S₂/N₂)*λ₂
end

function build_prob()
    u0 = [4.3e6*0.7 - 100.,100.,0.,0.,0.,0.,
            4.3e6*0.3 - 100.,100.,0.,0.,0.,0.,
            0.0,0.0]
    f = ODEFunction(twoSEStransmission;syms = [:S₁,:E₁,:I₁,:R₁,:W₁,:C₁,:S₂,:E₂,:I₂,:R₂,:W₂,:C₂,:C_sero₁,:C_sero₂])
    return ODEProblem(f,u0,(0,365))
end

"""
function simple_conv(incidence,kern)

Direct implementation of convolution [incidence ⋆ kern](t)
"""
function simple_conv(incidence,kern)
    z = similar(incidence)
    n = length(incidence)
    m = length(kern)
    for t = 1:n
        z[t] = 0.
        for s = 1:(t-1)
            if t-s <= m && t-s >= 1
                z[t] += kern[t-s]*incidence[s]
            end
        end
    end
    return z
end