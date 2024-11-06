function VFI2nd!(M::Model)
  # Note that many of these have an extra element in the age-dimension (shadow age). The policy in shadow ages are empty (since the dead have no choices)
  # while the shadow value function is just the altruistic continuation value
  gs = Vector{Pol}(undef,0)
  g = M.g
  val = M.V
  np = M.np
  mp = M.mp
  numdet = M.numdet
  Vknew = zeros(tuple(collect(size(val.Vk))..., np.viters+1))::Array{Float64,8}  # Final dimension is the value function iteration
  Vpnew = zeros(tuple(collect(size(val.Vp))..., np.viters+1))::Array{Float64,8}
  Vknew[:,:,:,:,:,:,:,1] = deepcopy(val.Vk) # Initial guess 
  Vpnew[:,:,:,:,:,:,:,1] = deepcopy(val.Vp)
  println("VFI!")

  # Preallocate interpllation arrays for 2nd stage problems:
  Vp2nd = [linear_interpolation((np.xphat_grd,),np.xphat_grd,extrapolation_bc=Throw()) for ixk = 1:np.nxk, ivk = 1:np.nv, iok = 1:np.no, iop = 1:np.no, is = 1:np.ns, ihp = 1:np.nh]
  tp2nd = [linear_interpolation((np.xphat_grd,),np.xphat_grd,extrapolation_bc=Throw()) for ixk = 1:np.nxk, ivk = 1:np.nv, iok = 1:np.no, iop = 1:np.no, is = 1:np.ns, ihp = 1:np.nh]
  
  # Preallocate interpolation arrays for problems
  tp_itp,bp_itp,discp = itp_gp_preallocate(np,g.gp);
  ck_itp,disck_itp,Vpnxt_itp = itp_gk_preallocate(np,g.gk,val.Vp);  
  
  ## Preallocate arrays for parents two stage problems:
  vtmp2nd = fill(-Inf64,np.nxphat,np.nh,np.ntp_choice2,np.nxk*np.nv*np.no*np.no*np.ns)
  vtmp2nd2 = fill(-Inf64,np.nxphat,np.nh,np.nxk*np.nv*np.no*np.no*np.ns)
  ttmp2nd2 = fill(-Inf64,np.nxphat,np.nh,np.nxk*np.nv*np.no*np.no*np.ns)
  vtmp1st = fill(-Inf64,np.ncp_choice2,np.nxk*np.nxp*np.nv*np.no*np.no*np.ns)
  
  iter = 0
  for _iter = numdet.previter+1:numdet.previter+np.viters
    iter += 1 # iter counts within each VFI, while _iter counts accross repeated VFI calls
    # Only used to store old policy function for distance calculation
    go = deepcopy(g)
    
    for iak = np.nyoung:-1:1      
      # Reset temporary value and policy functions  
      vtmp2nd .= -Inf64
      vtmp2nd2 .= -Inf64
      ttmp2nd2 .= -Inf64
      vtmp1st .= -Inf64
      
      g.gk.disc[:,:,:,:,:,:,iak,:] .= 0 

      # First solve kids problem at age iak, taking parents problem as given
      Vknxt_itp = itp_gp!(np,iak,g.gp,Vknew,Vpnew,tp_itp,bp_itp,discp,iter)
      Threads.@threads  for indx = 1:np.nxk*np.nbp*np.nh*np.nv*np.no*np.ns
        ixk,ibp,ihp,ivk,iok,is = Tuple(CartesianIndices((np.nxk,np.nbp,np.nh,np.nv,np.no,np.ns))[indx])
        Vknew[ixk,ibp,ihp,ivk,iok,is,iak,iter+1] = Valk!(np,mp,ixk,ibp,ihp,ivk,iok,is,iak,tp_itp,bp_itp,discp,Vknxt_itp,g,val)   
      end

      #  Then solve parents problem
      itp_gk!(np,iak,g.gk,Vpnew,ck_itp,disck_itp,Vpnxt_itp,iter);  
      # Solve 2nd stage problem (how much to give, for a given savings and housing choice:)
      @Threads.threads for indx = 1:np.nxk*np.nv*np.no*np.no*np.ns       
        ixk,ivk,iok,iop,is = Tuple(CartesianIndices((np.nxk,np.nv,np.no,np.no,np.ns))[indx])
        @views ValP2nd!(np,mp,ixk,ivk,iok,iop,is,iak,ck_itp,disck_itp,Vpnxt_itp,vtmp2nd[:,:,:,indx],vtmp2nd2[:,:,indx],ttmp2nd2[:,:,indx])          

        for ihp = 1:np.nh
          @views Vp2nd[ixk,ivk,iok,iop,is,ihp].itp.coefs .= vtmp2nd2[:,ihp,indx]
          @views tp2nd[ixk,ivk,iok,iop,is,ihp].itp.coefs .= ttmp2nd2[:,ihp,indx]
        end
      end

      # Solve first stage problem (consumption
      @Threads.threads for indx = 1:np.nxk*np.nxp*np.nv*np.no*np.no*np.ns
        ixk,ixp,ivk,iok,iop,is = Tuple(CartesianIndices((np.nxk,np.nxp,np.nv,np.no,np.no,np.ns))[indx])
        @views Vpnew[ixk,ixp,ivk,iok,iop,is,iak,iter+1] = ValP1st!(np,mp,ixk,ixp,ivk,iok,iop,is,iak,vtmp1st[:,indx],g,val,Vp2nd,tp2nd)
      end    
    end

    distks, distps, distances = find_VFI_distance(iter,Vknew,Vpnew,val)
    distk, distp = distks[1], distps[1] # Only take one lag for distance
    distance, idist = findmin(distances) # Find the smallest distance
       
    numdet.vdists[_iter] = distances[1] # Only store the distance for the first lag
    numdet.dist = distance
    numdet.previter = _iter
    
    VFI_display(g,go,Vknew,Vpnew,distance,distk,distp,iter)
    
    @views val.Vk = deepcopy(Vknew[:,:,:,:,:,:,:,iter+1])
    @views val.Vp = deepcopy(Vpnew[:,:,:,:,:,:,:,iter+1])
    
    push!(gs,deepcopy(g))
    if distance < M.np.vtol
      println("Converged")
      numdet.cnvrgd = true
      if idist == 1 # "Standard" convergence, one true fixed point
        numdet.cyclical = false
        numdet.cycle_length = idist
      else # Value functions converge cyclically (e.g., every other iteration has the same solution)
        numdet.cyclical = true
        numdet.cycle_length = idist
        @warn "Value functions converged to a cyclical fixed point. The cycle length is $idist"
      end

      break
    end
  end


  if numdet.cnvrgd == false
    throw("Value functions did not converge")
  end

  M.V=val
  M.g=g
  M.numdet = numdet
  
  return gs

end

function find_VFI_distance(iter,Vknew,Vpnew,val)
  distks = fill(0.0,iter)
  distps = fill(0.0,iter)  

  for i in 1:iter
    @views distks[i] = maximum(abs,(Vknew[:,:,:,:,:,:,:,iter+1] - Vknew[:,:,:,:,:,:,:,iter-i+1]))  
    @views distps[i] = maximum(abs,(Vpnew[:,:,:,:,:,:,:,iter+1] - Vpnew[:,:,:,:,:,:,:,iter-i+1]))  
  end
  distances = max.(distks,distps)
  return distks, distps, distances
end

function solvemodel2nd!(M::Model,np::NumPar,mp::ModPar)    
    
    @time   gs = VFI2nd!(M) # Solve decision problems  
    @time pan = simulate(M;gs) # Simulate households
    M.simpan = DataFrame(pan)  # Create dataframe of simulated panel
    @time M.moms = find_moments(np,mp,pan,np.switch_str["period"]) # Find aggregate moments

    ModelSolution_Display(M)
      

  end
  
function ModelSolution_Display(M)
  if M.numdet.cyclical == true
    @warn "Value functions converged to a cyclical fixed point. The cycle length is $(M.numdet.cycle_length)"
    display("The average moments of the simulated panels, by cycle, is:")
    vars = [:iok,:xk,:tp,:ivk]
    _df = combine(groupby(M.simpan,:cyc),vars.=>mean,renamecols=false)

    ## Calculate whether each indicator is within the 95% CI
    __df = M.simpan[M.simpan.cyc.==1,:]
    for var in vars
      stddev = [-1,1]*1.96*std(__df[!,var])./sqrt(nrow(__df))
      _df[!,string(var)*"_in95%CI"] = stddev[1] .<= (_df[!,var][1] .- _df[!,var]) .<= stddev[2] 
    end
    display(_df)
    display("The first set of columns is the means by cycle, while the second group shows whether each cycle mean is within the 95% CI of the first")
    display("Note that  difference between cycles include small sample difference, since these are housholds with different shock sequences as well")
    println("maximum percentage difference: ownership of kids: $(round(maximum(abs.((_df[!,:iok] .- _df[1,:iok])./_df[1,:iok]))*100;digits=3))%")
    println("maximum percentage difference: wealth    of kids: $(round(maximum(abs.((_df[!,:xk] .- _df[1,:xk])./_df[1,:xk]))*100;digits=3))%")
    println("maximum percentage difference: transfer from par: $(round(maximum(abs.((_df[!,:tp] .- _df[1,:tp])./_df[1,:tp]))*100;digits=3))%")
    
    println("")
  end

  display(M.moms[!,[:race,:medwealthyoung, :medwealthold, :owneryoung, :r2incyoung, :first_ownyoung, :wealthatpurchaseyoung, :transfrateyoung,:LTVatpurchaseyoung, :transferbuyersyoung]])
end
  

function findMb2nd(nstate,nchoice)
    mp, np = benchpar(nstate,nchoice)
    Mb = Model(mp,np)
    solvemodel2nd!(Mb,np,mp)

    return Mb,mp,np
end


function VFI_display(g,go,Vknew,Vpnew,distance,distk,distp,iter)
  # Calculate objects used to check for convergence (i.e., val function distance and policy function distances)
  dist_bk = go.gk.b- g.gk.b
  dist_dk = go.gk.disc - g.gk.disc
  dist_cp = go.gp.c - g.gp.c
  idist_cp = argmax(abs.(dist_cp))
  dist_tp = go.gp.t - g.gp.t
  idist_tp = argmax(abs.(dist_tp))
  dist_dp = go.gp.disc - g.gp.disc

  @views idist_vk = Tuple(argmax(abs.(Vknew[:,:,:,:,:,:,:,iter+1] - Vknew[:,:,:,:,:,:,:,iter])))
  @views idist_vp = Tuple(argmax(abs.(Vpnew[:,:,:,:,:,:,:,iter+1] - Vpnew[:,:,:,:,:,:,:,iter])))
  println(" ",(@sprintf "%3i" iter), " ", (@sprintf "%.10e" distance) ,
     " | k: ", (@sprintf "%.2e" distk), " @ ", join([@sprintf("%2i,", idist_vk[i]) for i in 1:7]),") " ,  
     " | p: ", (@sprintf "%.2e" distp), " @ ", join([@sprintf("%2i,", idist_vp[i]) for i in 1:7]),") " , 
     " | bk: ",  (@sprintf "%5.1f" maximum(abs,(dist_bk))), 
     " @ ("*join([@sprintf("%2i,", Tuple(argmax(abs.(dist_bk)))[i]) for i in 1:7]),"), " , 
     "dk: ",  (@sprintf "%1i" maximum(abs,(dist_dk))), 
     " | cp: ",  (@sprintf "%5.1f" dist_cp[idist_cp]),
     " @ ("*join([@sprintf("%2i,", Tuple(argmax(abs.(dist_cp)))[i]) for i in 1:8]),"), " , 
     "tp: ",  (@sprintf "%5.1f" dist_tp[idist_tp]), 
     " @ ("*join([@sprintf("%2i,", Tuple(idist_tp)[i]) for i in 1:8]),"), " , 
     "dp: ",  (@sprintf "%1i" maximum(abs,(dist_dp)))
     )
end
