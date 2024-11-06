function isowner(ih,np::NumPar)
  @assert ih>0
  @assert ih<=np.nh

  owner = true
  if ih <= np.nhr
    owner = false
  end
  return owner

end

function BorrCons(h,own::Bool,price,ia::Int,mp::ModPar)
  LTV = - mp.LTV[ia]*price*h*own

  BorrCons = LTV
end


function hexp_func(h::Float64,own::Bool,ia::Int64,io::Int64,price::Float64,np::NumPar,mp::ModPar)
  if own == false
    hexp = mp.κ*price*h
  else
    hexp = price*h
  end

  if own == true && io <= np.nhr # Buyers must pay extra..
    hexp = hexp + mp.mb[ia]*price*h
  elseif own == false && io > np.nhr # Sellers must pay
    hexp = hexp + mp.ms[ia]*price*np.h_grd[io]
  end

  return hexp
end

function BCp_netshock(x::Float64,c::Float64,t::Float64,h::Float64,own::Bool,ia::Int64,io::Int64,price::Float64,np::NumPar,mp::ModPar)
  @assert ia>np.nyoung
  hexp = hexp_func(h,own,ia,io,price,np,mp)

  b = x - c - t - hexp + wp_netshock(ia,np)
  
  return b
end

function BCk(x::Float64,c::Float64,t::Float64,h::Float64,own::Bool,ia::Int64,io::Int64,price::Float64,np::NumPar,mp::ModPar,iinc)
  @assert ia≤np.nyoung

  hexp = hexp_func(h,own,ia,io,price,np,mp)
  b = x - c - t - hexp + wk(ia,iinc,np)
    
  return b
end

function xk′_scalar(b,h,own::Bool,iak,np::NumPar,mp::ModPar,is,iδ)
    @assert iak < np.nyoung && iak > 0
    xkn = b*(1+rr(b,iak,h,own,np.p_grd[is],mp)) + h*own*(1-mp.δ-np.δ_grd[iδ])*np.p_grd[is]
end

function xkdead′_scalar(b,h,own::Bool,iak,np::NumPar,mp::ModPar,is,iδ,iν) # In the last period, next period wealth of the current kid depends on his transitory income shock as a parent.
  @assert iak == np.nyoung 
  xkn = b*(1+rr(b,iak,h,own,np.p_grd[is],mp)) + h*own*(1-mp.δ-np.δ_grd[iδ])*np.p_grd[is] + wp_shock(iak+1,iν,np)
end

function xp′_scalar(b,h,own::Bool,iap,np::NumPar,mp::ModPar,is,iδ,iν)
    @assert iap > np.nyoung && iap < np.na
    xpn = b*(1+rr(b,iap,h,own,np.p_grd[is],mp)) + h*own*(1-mp.δ-np.δ_grd[iδ])*np.p_grd[is] + wp_shock(iap+1,iν,np)
end

function xpdead′_scalar(b,h,own::Bool,iap,np::NumPar,mp::ModPar,is,iδ)
  @assert iap ==  np.na
  xpn = b*(1+rr(b,iap,h,own,np.p_grd[is],mp)) + (h*own*(1-mp.δ-np.δ_grd[iδ])*np.p_grd[is])*(1-mp.ms[iap])
end


function rr(bk,ia::Int,h,own::Bool,price,mp::ModPar)
  if bk>.0
    r = mp.rf
  else
    r = mp.rf + mp.rm[ia]
    if bk < - mp.PMIlim*price*h*own 
      r += mp.PMIprem
    end
    
  end

  return r

end

function wk(iak::Int64,ivk::Int64,np::NumPar)
  @assert iak <= np.nyoung
  wk = np.y_grd[iak]*np.v_grd[iak,ivk]
end

# Parents wage
function wp(iap::Int64,iνp::Int64,np::NumPar) 
  @assert iap > np.nyoung
  wk = np.y_grd[iap]*np.ν_grd[iap,iνp]
end

# The Income shock from parents income shock (only used in value function iteration, so that parents income isn't a state)
function wp_shock(iap::Int64,iνp::Int64,np::NumPar)
  @assert iap > np.nyoung
    wk = np.y_grd[iap]*(np.ν_grd[iap,iνp]-1.0)
end

# Parents wage, net of shock (only used in value function iteration)
function wp_netshock(iap::Int64,np::NumPar)
  @assert iap > np.nyoung
  wk = np.y_grd[iap]
end

function util(c::Float64,h::Float64,own::Bool,mp)
  @assert c > 0.0 && h>0.0
  if own == true
    ch = c*(h*(mp.χ)/c )^(mp.ξ) # Note how you've divied by c (c^(1-ξ)*h^ξ = c*(h/c)^ξ) to get rid of one exponent ^ calculation
  else
    ch = c*(h/c)^(mp.ξ)
  end
  util = (ch^(1.0-mp.γ) - 1.0)/(1.0 - mp.γ)

  return util::Float64
end