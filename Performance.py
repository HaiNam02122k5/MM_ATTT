import NumberTheory
import cProfile

# Profile func1
profiler = cProfile.Profile()
profiler.enable()
NumberTheory.primitiveRoot(117809)
profiler.disable()
profiler.print_stats(sort='time')  # In kết quả, sort theo thời gian
