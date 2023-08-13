
# def myMap(*args, **kwargs):
#   return map(*args, **kwargs)


def myMap(*args, **kwargs):
  """
  Applies args[0] to each element of args[1]. Returns a list of the results.
  :param args:
  :param kwargs:
  :return:
  """
  return [args[0](x) for x in args[1]]


# from joblib import Parallel, delayed
# def myMap(*args, **kwargs):
#   return Parallel(n_jobs=5, backend="loky")(delayed(args[0])(arg) for arg in args[1])