include("../HomeownershipBankMomDad.jl")


#######################
## This script runs all Julia code to reproduce the entire paper, including appendix.
## Note that the structural estimation is slow. For reference, it takes about one week on my computational cluster using 30 threads.
#######################
## Set "parameters"
pars = [:η, :χ, :ho] # Estimated parameters
targ_moms =  [:owneryoung, :r2incyoung, :wealthatpurchaseyoung,:transfrateyoung]
nontarg_moms = [:medwealthyoung, :medwealthold, :medwealthpgradyoung, :first_ownyoung, :ownerall, :mortgyoung,:LTVatpurchaseyoung, :tp2wealthpyoung, :transferbuyersyoung,:transfraterenteryoung,:transfrateowneryoung]
nstate = 65 # Grid points in state
nchoice = 145 # Grid points in choice
Nest = 1500 # How many iterations in each step of the estimation
Nshrinks = 2 # How many times to "shrink" the estimation area (if 0, you don't do it at all)
store = true # true to store results to .pdf/.tex



#######################
## Structural Estimation
#######################

for shrinkstep = 0:Nshrinks # Shrink the search space iteralatively
    if shrinkstep == 0 # Initialize search space
        par_lims = par_lims_ini()
    else # Shrink search space, based on the last Nest iterations
        ilo = 1 + Nest*(shrinkstep-1)
        ihi = Nest*shrinkstep
        par_lims = shrink_searchspace(ilo,ihi,nstate,nchoice,par_lims);
    end
    GlobalSearch(par_lims, nstate,nchoice,targ_moms;nsim = (Nest*(shrinkstep+1)),start=(1+Nest*(shrinkstep)))
end

## Plot estimation
df_boot, globest,dat = LoadEstimation(nstate,nchoice,targ_moms,pars);
PlotGlobalEstimation(dat,globest,pars,targ_moms;store=true)

#######################
## Model fit etc
#######################
# (re-)run the model with the best parameters so that you have the estimated panel
Mb,mp,np = findMb2nd(nstate,nchoice); # Need to re-run model to obtain simulated panel
leftjoin!(Mb.moms,DataFrame(type=["Data","Model"],smm_obj=[NaN,smm_obj(Mb.moms[Mb.moms.type .== "Model",:],Mb.moms[Mb.moms.type .== "Data",:],targ_moms)]),on=:type) # Calculate the SMM objective function and store it in the Mb.moms dataframe

# Make sure that model output is the same as the one from the estimation
for par in [:η, :χ, :ho] # Ensure that the parameters are (almost) the same as the best from the estimation
    @assert abs(getfield(mp,par) - sort!(globest,:smm_obj)[!,par][1] ) .< 1e-5  "$par is not the same in 'mp' and 'est'"
end
for mom in targ_moms # Check that the data are exactly the same as they should be!
    @assert Mb.moms[Mb.moms.type.=="Data",mom] == dat[!, mom] "Data and model do not match"
end
@assert abs(sort!(globest.smm_obj)[1]) == Mb.moms[Mb.moms.type .== "Model",:smm_obj][1]

output_esttable(df_boot,pars,targ_moms,Mb,store=store)
output_extramoments(Mb,nontarg_moms,store=store)
pChetty = ChettyEventStudy(Mb,0.75,store=store)
#######################
## Post-Estimation Analaysis - Impact of Altruism on Homeownership
#######################
## Contribution of Altruism to Homeownership
Mn = altrcontrhousing(Mb.mp,Mb.np);
nontarg_moms_notransf = [x for x in nontarg_moms if x ∉ [:tp2wealthpyoung, :transferbuyersyoung,:transfraterenteryoung,:transfrateowneryoung]] # No reason to report transfer moments in this table
output_noalt(Mb,Mn,targ_moms,nontarg_moms_notransf;store=store)

## Policy function plts
plot_policyfuncs(Mb,Mn,store=store)

## Endogenous Prices
# First: Solve model for many prices
Mb,mp,np = findMb2nd(nstate,nchoice)
Hbar = Mb.moms[Mb.moms.type .== "Model",:].ownerall[1]
Pbench = mp.price
pgrd = range(start=Mb.mp.price*0.8,stop=Mb.mp.price,length=10)

dfdemand = SolveModelDifferentPrices(Mb,pgrd) # Solve model without altruism for various prices
##  Calculate the supply for each price and supply elasticity
α1s = [1.0, 3.0, Inf]
α0s = [Hsupply_intercept(Hbar,Pbench,α1) for α1 in α1s]

supplymat = fill(0.0,length(pgrd),length(α1s),) # Rows are alphas, cols are prices
for iα in eachindex(α1s)
    α1 = α1s[iα]
    α0 = α0s[iα]
    for (ip,price) in enumerate(pgrd)
        supplymat[ip,iα] = supply_func(α0,α1,price)
    end
end

# Find equilibrium (minimizes distance between supply and demand)
eqlbm = DataFrame()
for iα in eachindex(α1s)
    ip = argmin(abs.(dfdemand.demand - supplymat[:,iα]))
    supply=supplymat[ip,iα]
    if α1s[iα]==Inf
        ip = argmin(abs.(pgrd .- mp.price))
        supply = dfdemand.demand[ip] # If perfectly elastic, supply is whatever demand isempty
    end
    append!(eqlbm,DataFrame(α1=α1s[iα],α0=α0s[iα],supply=supply,demand=dfdemand.demand[ip],price=dfdemand.price[ip],ip=ip))
end

sort!(eqlbm,:α1)
@assert eqlbm.ip[1] > 1 "Need eqlbm to be on the inside of the price grid: Decrease the lowest price"
@assert maximum(eqlbm.ip[1:end-1]) < length(pgrd) "Need eqlbm to be on the inside of the price grid: Top of the grid is too sparse"
@assert length(unique(eqlbm.ip))==length(α1s) "Check that you get different prices each supply: The grid is too sparse"
#
plot(dfdemand.price,dfdemand.demand,label="Demand w/o Altruism",xlabel="Price",ylabel="Housing")
for (iα,α1) in enumerate(α1s)
    
    if α1 == Inf
        vline!([mp.price],label="Supply (α1=$α1)",marker=:circle)
    else 
        plot!(dfdemand.price,supplymat[:,iα],label="Supply (α1=$α1)",marker=:circle,)
    end
end
display(plot!())

aggmoms = [ :α1,:price,:ownerall]
dfp = innerjoin(dfdemand,eqlbm,on=:price,makeunique=true)

output_endoprice(dfp,Mb,targ_moms,aggmoms,α1s;store=true)




#######################
## Post-Estimation Analaysis - Policy, Illiquidity, Parental Wealth
#######################
Mb,mp,np = findMb2nd(nstate,nchoice)
out = policyeffects(Mb)
bypwealth_moms = [:owneryoung_prich, :owneryoung_pmiddle,:owneryoung_ppoor]
output_policy(out,targ_moms,bypwealth_moms;store=store)

## MPC analysis
create_MPC_table(Mb,Mn,store)

## Adjustment costs
Mb,mp,np = findMb2nd(nstate,nchoice)
Mn = solve_adjustmentcost(Mb)
df1,df2,dfall = quant_adjustmentcost(Mb,Mn)
adjustmentcost_table(df2;store=true)
plot(adjustmentcost_plot(dfall,store=store),size=(500,500))

## Black-White
Mb,mp,np = findMb2nd(nstate,nchoice)
BW = blackwhitegap(Mb)
output_BW(BW,targ_moms,[];store=store)

#############################
## Appendix
#############################
Mb,mp,np = findMb2nd(nstate,nchoice)
Mn = altrcontrhousing(mp,np)

Mbν,Mnν,Mbs,Mns = robustness(nstate,nchoice)
output_robust(Mb,Mn,Mbν,Mnν,Mbs,Mns,targ_moms,nontarg_moms;store=store)

## Create latex parameter values.
open("tabfig/parameters.tex", "w") do f
    write(f,"\\newcommand{\\parxi}{$(mp.ξ)}\n")
    write(f,"\\newcommand{\\parbeta}{$(mp.β)}\n")
    write(f,"\\newcommand{\\pargamma}{$(mp.γ)}\n")
    write(f,"\\newcommand{\\parq}{$(mp.κ)}\n")
    write(f,"\\newcommand{\\pardelta}{$(mp.δ)}\n")
    write(f,"\\newcommand{\\parrf}{$(mp.rf)}\n")
    write(f,"\\newcommand{\\parrm}{$(mp.rm[1])}\n")
    write(f,"\\newcommand{\\parrpmi}{$(mp.PMIprem)}\n")
    write(f,"\\newcommand{\\parPMIlim}{$(mp.PMIlim)}\n")
    write(f,"\\newcommand{\\parLTV}{$(mp.LTV[1])}\n")
    write(f,"\\newcommand{\\parms}{$(mp.ms[1])}\n")
    write(f,"\\newcommand{\\parmsperc}{$(round(Int,mp.ms[1]*100))}\n")
    write(f,"\\newcommand{\\parmb}{$(mp.mb[1])}\n")
    write(f,"\\newcommand{\\parmbperc}{$(round(Int,mp.mb[1]*100))}\n")
    write(f,"\\newcommand{\\parprice}{$(round(mp.price;digits=2))\$\\times h_o\$}\n")
    write(f,"\\newcommand{\\parhr}{$(np.h_grd[1])}\n")    

    write(f,"\\newcommand{\\parNdyn}{$(np.Ndyn)}\n")    
    write(f,"\\newcommand{\\parNgen}{$(np.Ngen)}\n")    
    write(f,"\\newcommand{\\parnstate}{$(nstate)}\n")    
    write(f,"\\newcommand{\\parnchoice}{$(nchoice)}\n")    
    write(f,"\\newcommand{\\parNest}{$(Nest)}\n")    
    write(f,"\\newcommand{\\parNshrinks}{$(Nshrinks)}\n")    
    

end
