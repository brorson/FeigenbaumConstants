# This program computes Feigenbaum's alpha constant using
# Feigenbaum's functional equation derived via renormalization
# group methods.  This version uses a Taylor's expansion for
# the function g(z).  The expansion is evaluated on a set of
# Chebyshev nodes.
#
# Started: Nov 2017 -- Stuart Brorson, sdb@cloud9.net

using ForwardDiff
using JLD

#++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Helper fcns
    

#-------------------------------------------------------
# Now create the Chebyshev points
function generate_cheb_points(N)
  #println("---->  Entering generate_cheb_points\n")
  local k = big(1.0):big(N)
  local bigpi = big(pi);
  zi = cos.(bigpi*(2*k-1)./(4*N));
  #println(zi)
  #println("<----  Leaving generate_cheb_points\n")
  return zi;
end


#-------------------------------------------------------
# Define universal fcn g(z, a).  It is expressed as a
# Taylor's series, coefficients a[i].  Only even powers
# are needed.  Use Horner's rule to perform sum.
function g{T<:Real}(z, a::Vector{T})
  #println("  ---->  Entering g, z = $z\n") 
  local s = zero(z);
  for i=length(a):-1:1
    s = z*z*(s + a[i]);
  end
  #println("  <----  Leaving g, s = $s\n") 
  return (one(z)+s)
end

#---------------------------------------------------------
# Define relation obeyed by g(z).  Name it f(z, a)
function f{T<:Real}(z, a::Vector{T})
  #println("---->  Entering f, a = $a\n")
  local myone = one(z);
  local alpha = myone/g(myone, a);
  r = g(z, a) - alpha*g(g(z/alpha, a), a);
  #println("<----  Leaving f\n")
  return r
end


#--------------------------------------------------------
# Define fcn which takes gradient of f w.r.t. a.
# Returns row vector which is gradient of f w.r.t. a,
# evaluated at z, a.
function gradf(z, a::Vector{BigFloat})
  # println("---> Entered gradf\n")
  # Specialize to position z.  Result is fcn of coeffs a only.
  fbeta(a) = f(z, a);

  # Compute gradient w.r.t. beta and return it
  gf1 = ForwardDiff.gradient(fbeta, a)
  # println("<--- Leaving gradf\n")
  return gf1
end


# ------------------------------------------------------------
# The main event -- a computation of alpha using Newton's method.
function compute_alpha(Numdigs, Npts, betain)

  # First set up variables, including BigFloat settings
  # Setting precision to 5x Numdigs seems to keep Newton's method
  # from wandering away.  Smaller values don't.  More fidding might
  # improve things here.
  setprecision(3*Numdigs);
  
  # Set stopping tol for Newton's method so I get the number of digits
  # I want.
  const tol = BigFloat(10.0)^-(Npts);

  zi = generate_cheb_points(Npts);      # Grid of sample points zi

  fn = zeros(BigFloat, Npts);             # Function f vector
  Jn = zeros(BigFloat, Npts, Npts);        # Jacobian
 
  # Need expansion coeffs beta.  Start with the one passed in, betain.  
  # N should always be larger than length(betain), so we are adding to the 
  # betas we computed last time.
  betan = zeros(BigFloat, Npts);
  for i = 1:length(betain)
    betan[i] = betain[i];
  end

  # println("-----------------------------\n")
  sttim = now();  # Keep track of loop timing.

  # Now enter Newton loop.  
  for cnt = 1:10

    # Compute new f and Jacobian upon each iteration
    #println("Computing fn and Jn...")
    for i = 1:Npts
      fn[i] = f(zi[i], betan);
      Jn[i,:] = transpose(gradf(zi[i], betan));
    end

    #println(fn)
    #println(Jn)

    # Now compute step to take, then take it.
    #println("Computing step sn...")
    sn = Jn\fn;
    #println("Step = $sn\n")
    betanp1 = betan - sn;
    #println("betanp1 = $betanp1\n")

    nowtim = now();  # Keep track of loop timing
    println("After Newton iteration = $cnt, time elapsed = $(nowtim-sttim)\n")

    # Write out betanp1 and Npts values in case we need to restart the
    # calculation
    save("./betan.jld", "betan", betanp1, "Npts", Npts)

    # Check for convergence
    if (norm(sn) < tol)
      println("\nDone!  Converged after $cnt iterations.\n")
      # compute alpha
      myalpha = big(1.0)/g(big(1.0), betanp1);
      return myalpha, betanp1;
    end

    # Move values back
    betan = betanp1

  end
  println("!!!!!!!  Didn't converge!  !!!!!!!!\n")

end


#--------------------------------------------------
function iterate_alpha()
  # This is a runner -- call this fcn.  This fcn calls
  # compute_alpha asking for different number of digits 
  # at different precision/tolerance levels.
  # It calls with the beta coefficients computed in the 
  # last iteration as the starting point to help Newton's 
  # method from wandering away.

  # Seed betan, either with known start values, or with
  # the values computed at the last Newton iteration.
  println("Use betan value from file [y/N]? ==>")
  resp = readline(STDIN)
  if (resp == "")
    betan = zeros(BigFloat, 3);
    betan[1] = big(-1.5);
    betan[2] = big(1e-1);
    betan[3] = big(2e-2);
    N0 = 4;
  else
    betan = load("./betan.jld", "betan")
    N0 = load("./betan.jld", "Npts")
  end
 
  # Step up number of pts to compute on.  Change this to get
  # more digits
  for Npts = N0:2*N0:2500
    println("=================================\n")
    println("Npts = $Npts\n")
    Numdigs = Int(floor(2*Npts));

    # Call the function which computes alpha at this level
    # of precision.
    myalpha, betan = compute_alpha(Numdigs, Npts, betan)
    
    # Print out myalpha
    a = string(myalpha);
    println("myalpha = $(a)\n")

    f = open("alpha.dat", "a");
    write(f, "==========\nNpts = $Npts\n$(a)\n")
    close(f)

    # Compare against the version from 
    # http://www.plouffe.fr/simon/constants/feigenbaum.txt
    oeisalpha = oeis_alpha(3*Numdigs);
    diff = Float64(abs(myalpha) - abs(oeisalpha));
    println("diff = $diff\n")
  end

end

#------------------------------------------------
function oeis_alpha(Numdigs)
  # Return up to 1018 digits of alpha from
  # http://www.plouffe.fr/simon/constants/feigenbaum.txt

  if (Numdigs <= 1018)
    setprecision(Numdigs);
  else
    setprecision(1018);
  end

  str = "-2."
  str = str*"5029078750958928222839028732182157863812713767271499773361920567"
  str = str*"7923546317959020670329964974643383412959523186999585472394218237"
  str = str*"7785445179272863314993372578112163594879503744781260997380598671"
  str = str*"2397117373289276654044010306698313834600094139322364490657889951"
  str = str*"2205843172507873377463087853424285351988587500042358246918740820"
  str = str*"4281700901714823051821621619413199856066129382742649709844084470"
  str = str*"1008054549677936760888126446406885181552709324007542506497157047"
  str = str*"0475419932831783645332562415378693957125097066387979492654623137"
  str = str*"6745918909813116752434221110130913127837160951158341230841503716"
  str = str*"4997020224681219644081216686527458043026245782561067150138521821"
  str = str*"6449532543349873487413352795815351016583605455763513276501810781"
  str = str*"1948369459574850237398235452625632779475397269902012891516645793"
  str = str*"9420198920248803394051699686551494477396533876979741232354061781"
  str = str*"9896112494095990353128997733611849847377946108428833293833903950"
  str = str*"9008914086351525626803381414669279913310743349705143545201344643"
  str = str*"4264752001621384610729922641994332772918977769053802596851"

  alpha = parse(BigFloat, str)

end
