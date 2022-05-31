
# def myMap(*args, **kwargs):
#   return map(*args, **kwargs)


def myMap(*args, **kwargs):
  return [args[0](x) for x in args[1]]


# from joblib import Parallel, delayed
# def myMap(*args, **kwargs):
#   return Parallel(n_jobs=5, backend="loky")(delayed(args[0])(arg) for arg in args[1])