##################
## kids decision problems
#################
function evalk(np::NumPar,mp::ModPar,ck::Float64,ihk::Int64,ixk,ibp,ihp,ivk,iok,is,iak,tp_itp,bp_itp,discp,Vnxt)
  price = np.p_grd[is]

  hk = np.h_grd[ihk]

  ownk = isowner(ihk,np)
  bk = BCk(np.xk_grd[ixk],ck,0.0,hk,ownk,iak,iok,price,np,mp,ivk)
  bkn = bk*(1.0 + rr(bk,iak,hk,ownk,price,mp))
  ownp = isowner(ihp,np)

  if ck > 0.0 && bk >= BorrCons(hk,ownk,price,iak,mp)
    contval = util(ck,hk,ownk,mp)

    if iak < np.nyoung
      for isn = 1:np.ns
        if np.Πs[is,isn] > 0
        for iνn = 1:np.nν
            for iδp = 1:np.nδ
            xpn = xp′_scalar(np.bp_grd[ibp],np.h_grd[ihp],ownp,iak+np.nyoung,np,mp,isn,iδp,iνn)
              for iδk = 1:np.nδ
                xkn = xk′_scalar(0.0,hk,ownk,iak,np,mp,isn,iδk) + bkn
                for ivn=1:np.nv
                # Loop over possible next-period housing choices of the parent
                for ihpn = 1:np.nh
                  prob = discp[ivn,ihk,ihp,isn,ihpn](xkn,xpn)::Float64
                  if prob > 0.0
                    tpn = (1.0-mp.τt)*tp_itp[ivn,ihk,ihp,isn,ihpn](xkn,xpn)::Float64
                    bpn = bp_itp[ivn,ihk,ihp,isn,ihpn](xkn,xpn)::Float64
                    contval += mp.β*np.Πs[is,isn]*np.Πv[iak,ivk,ivn]*np.Πν[iak+np.nyoung+1,iνn]*np.Πδ[iδk]*np.Πδ[iδp]*prob*Vnxt[ihpn,ivn,ihk,isn](xkn+tpn,bpn)::Float64
                  end
                end
              end
            end
          end
          end
        end
      end
    else
      for ixn = 1:np.Θnxk
        for ivn = 1:np.nv
          for iνn = 1:np.nν
            for isn = 1:np.ns
              for iδp = 1:np.nδ
                xpn = xpdead′_scalar(np.bp_grd[ibp],np.h_grd[ihp],ownp,iak+np.nyoung,np,mp,isn,iδp)*(1.0 - mp.τb)
                for iδk = 1:np.nδ
                  xkn = xkdead′_scalar(0.0,hk,ownk,iak,np,mp,isn,iδk,iνn) + bkn
                  contval += mp.β*np.Πs[is,isn]*np.ΠΘ[ixn,ivn,ixk,ivk]*np.Πν[iak+1,iνn]*np.Πδ[iδk]*np.Πδ[iδp]*Vnxt[mp.ihini,ivn,ihk,isn](np.θx_grd[ixn,ivn],xkn + xpn)::Float64
                end
              end
            end
          end
        end
      end
    end

  else
    contval = -1.0e40
  end

  return contval#, ck
end

function Valk!(np::NumPar,mp::ModPar,ixk,ibp,ihp,ivk,iok,is,iak,tp_itp,bp_itp,discp,Vnxt_itp,
  g,val)
  
  xk = np.xk_grd[ixk]
  price = np.p_grd[is]

  _val = -Inf
  _itmp = 0

  for ihk in 1:np.nh
      hk = np.h_grd[ihk]
      ownk = isowner(ihk,np)

      valmax = - Inf
      imax = 1
      
      for ic in 1:np.nck_choice
          ck  = np.ck_grd[ic]
          _tmp = evalk(np,mp,ck,ihk,ixk,ibp,ihp,ivk,iok,is,iak,tp_itp,bp_itp,discp,Vnxt_itp)
          if _tmp > valmax # Update best value and policy
            valmax = _tmp
            imax = ic
          end
          if _tmp ≈ -1.0e40 # If no feasible choices remain, then exit the loop..!
            break
          end
      end

      ## Update conditional (on housing) policy and value functions 

      val.Vk_disc[ixk,ibp,ihp,ivk,iok,is,iak,ihk] = valmax
      g.gk.c[ixk,ibp,ihp,ivk,iok,is,iak,ihk] = ck = np.ck_grd[imax]
      g.gk.b[ixk,ibp,ihp,ivk,iok,is,iak,ihk] = BCk(xk,ck,0.0,hk,ownk,iak,iok,price,np,mp,ivk)
      g.gk.h[ixk,ibp,ihp,ivk,iok,is,iak,ihk] = hk      

      if valmax>_val # stores the best value discrete value
        _val = valmax
        _itmp = ihk
      end        
  end

  # Find the best discrete choice
  for ihk = 1:np.nh 
    if _itmp == ihk
      g.gk.disc[ixk,ibp,ihp,ivk,iok,is,iak,ihk] = 1
    end
  end

  return _val

end


##################
## Parents
##################
function ValP2nd!(np::NumPar,mp::ModPar,ixk,ivk,iok,iop,is,iak,gck_itp,gdisck_itp,Vnxt_itp,
  _vtmp2nd,_vtmp2nd2,_ttmp2nd2)
  price = np.p_grd[is]
  for ihp = 1:np.nh
    hp = np.h_grd[ihp]
    for (ixp,xphat) in enumerate(np.xphat_grd)
      for (it,tp) in enumerate(np.tp_grd2)
        _vtmp2nd[ixp,ihp,it] = evalp2nd(np,mp,xphat,tp,ihp,hp,ixk,ivk,iok,iop,is,iak,price,gck_itp,gdisck_itp,Vnxt_itp)
      end
      _vtmp2nd2[ixp,ihp],itmp = findmax(@view _vtmp2nd[ixp,ihp,:]) 
      _ttmp2nd2[ixp,ihp] = np.tp_grd2[itmp]
    end
  end  
end


function ValP1st!(np::NumPar,mp::ModPar,ixk,ixp,ivk,iok,iop,is,iak,vtmp1st,g,val,Vp2nd,tp2nd)
  xp = np.xp_grd[ixp]
  price = np.p_grd[is]
  for ihp = 1:2
    hp = np.h_grd[ihp]; ownp = isowner(ihp,np)
    for (icp,_cp) in enumerate(np.cp_grd2)
        xphat =  xp-_cp
        if _cp > 0 && xphat>np.xphat_grd[1] 
            vtmp1st[icp] = util(_cp,hp,ownp,mp) + Vp2nd[ixk,ivk,iok,iop,is,ihp](xphat)
        end    
    end
    v1,icp = findmax(vtmp1st)
    g.gp.c[ixk,ixp,ivk,iok,iop,is,iak,ihp] = cpstar = np.cp_grd2[icp]          
    g.gp.t[ixk,ixp,ivk,iok,iop,is,iak,ihp] = tpstar = tp2nd[ixk,ivk,iok,iop,is,ihp](xp-cpstar)
    g.gp.b[ixk,ixp,ivk,iok,iop,is,iak,ihp] = BCp_netshock(xp,cpstar,tpstar,hp,ownp,iak+np.nyoung,iop,price,np,mp)
    g.gp.h[ixk,ixp,ivk,iok,iop,is,iak,ihp] = np.h_grd[ihp]
    val.Vp_disc[ixk,ixp,ivk,iok,iop,is,iak,ihp] = v1
  end

  @views _val, imax_disc = findmax(val.Vp_disc[ixk,ixp,ivk,iok,iop,is,iak,:])
  for ihp = 1:np.nh # Update policies
    if imax_disc == ihp
      g.gp.disc[ixk,ixp,ivk,iok,iop,is,iak,ihp] = 1
    else
      g.gp.disc[ixk,ixp,ivk,iok,iop,is,iak,ihp] = 0
    end
  end
  
  return _val
end

function evalp2nd(np::NumPar,mp::ModPar,xphat::Float64,tp::Float64,ihp::Int64,hp,ixk,ivk,iok,iop,is,iak,price,ck_itp,disck,Vnxt)
  iap = iak + np.nyoung
  xp = xphat
  xk = np.xk_grd[ixk]
  cp = 0.0 # Placeholder for finding bp :)  
  xkcash = xk + tp*(1.0-mp.τt)

  ownp = isowner(ihp,np)
  bp = BCp_netshock(xp,cp,tp,hp,ownp,iap,iop,price,np,mp) 
  bpn = bp*(1.0 + rr(bp,iap,hp,ownp,price,mp))

  
  if (iap == 30 && tp > 0 && np.FinalTransfer == false )
    bp = -Inf
  end


  if bp >= BorrCons(hp,ownp,price,iap,mp) 
    contval = 0.0 

    for ihk = 1:np.nh
      prob = disck[ihp,ivk,iok,is,ihk](xkcash,bp)
      if prob > 0

        hk = np.h_grd[ihk]
        ownk = isowner(ihk,np)
        ck = ck_itp[ihp,ivk,iok,is,ihk](xkcash,bp)
        bk = BCk(xkcash,ck,0.0,hk,ownk,iak,iok,price,np,mp,ivk)
        bkn = bk*(1.0 + rr(bk,iak,hk,ownk,price,mp))

        if bk < BorrCons(hk,ownk,price,iak,mp)
          ck = ck - (bk -BorrCons(hk,ownk,price,iak,mp))
          bk = BCk(xkcash,ck,0.0,hk,ownk,iak,iok,price,np,mp,ivk)
        end
        contval += mp.η*util(ck,hk,ownk,mp)*prob
        if iak < np.nyoung
          for isn = 1:np.ns
            if np.Πs[is,isn] > 0
                for iν = 1:np.nν
                  for iδp = 1:np.nδ
                    xpn = bpn + xp′_scalar(0.0,hp,isowner(ihp,np),iap,np,mp,isn,iδp,iν)
                    for iδk = 1:np.nδ
                      xkn = bkn + xk′_scalar(0.0,hk,isowner(ihk,np),iak,np,mp,isn,iδk)                          
                      for ivn = 1:np.nv
                      contval += mp.β*np.Πs[is,isn]*np.Πv[iak,ivk,ivn]*np.Πν[iap+1,iν]*np.Πδ[iδk]*np.Πδ[iδp]*prob*Vnxt[ivn,ihk,ihp,isn](xkn,xpn)
                    end
                  end
                end
              end
            end
          end
        else
          for isn = 1:np.ns
            if np.Πs[is,isn] > 0.0
              for iν = 1:np.nν
                for iδp = 1:np.nδ
                  xpn = (bpn + xpdead′_scalar(0.0,hp,isowner(ihp,np),iap,np,mp,isn,iδp))*(1.0 - mp.τb)
                  for iδk = 1:np.nδ
                    for ixn = 1:np.Θnxk
                      xkn = bkn + xkdead′_scalar(0.0,hk,isowner(ihk,np),iak,np,mp,isn,iδk,iν)
                      for ivn = 1:np.nv                   
                        contval += mp.η*mp.β*np.Πs[is,isn]*np.Πν[iak+1,iν]*np.ΠΘ[ixn,ivn,ixk,ivk]*np.Πδ[iδk]*np.Πδ[iδp]*prob*Vnxt[ivn,mp.ihini,ihk,isn](np.θx_grd[ixn,ivn],xkn + xpn)
                      end
                    end
                  end
                end
              end
            end
          end
        end
      end
    end
  else
    contval = -1.0e40
  end

  return contval
end
