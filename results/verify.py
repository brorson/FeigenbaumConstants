

def verify():

  fname_test = raw_input("Enter the file under test (default = alpha_me.txt): ")
  if (fname_test == ""):
    fname_test = "alpha_me.txt"
  f_test = open(fname_test, "r")

  fname_ref = raw_input("Enter the reference file (default = alpha_oeis.txt): ")
  if (fname_ref == ""):
    fname_ref = "alpha_oeis.txt"
  f_ref = open(fname_ref, "r")

  test_obj = f_test.readlines()
  test = test_obj[0]
  ref_obj = f_ref.readlines()
  ref = ref_obj[0]

  # Count matches
  Ntest = len(test)
  Nref = len(ref)
  N = min(Ntest, Nref)
  cnt = 0
  s = ""
  for i in range(N):
    if (test[i] == ref[i]):
      s = s+test[i]
      cnt = cnt+1
    else:
      break
  
  print(fname_test+" agrees with "+fname_ref+" to "+str(cnt)+" digits.")
  print(fname_test+" has "+str(Ntest)+" digits.")
  print(fname_ref+" has "+str(Nref)+" digits.")
  print("Digits of agreement: alpha = \n"+s)

if __name__ == "__main__":
  verify()
