def show(iteration, total, prefix = '', suffix = '', decimals = 3, 
         length = 50, fill = '█', print_end = "\r"):
    """
    Call in a loop to create terminal progress bar. Taken from [1].    

    Parameters
    ----------
    iteration : int
        Current iteration.
    total : int
        Total iteration.
    prefix : str, optional
        Prefix string. The default is ''.
    suffix : str, optional
        Suffix string. The default is ''.
    decimals : int, optional
        Positive number of decimals in percent complete. The default is 1.
    length : int, optional
        Character length of bar. The default is 20.
    fill : str, optional
        Bar fill character. The default is '█'.
    print_end : str, optional
        End character. The default is "\r".

    Returns
    -------
    None.
    
    References
    -------
    .. [1] https://stackoverflow.com/a/34325723

    """
    
    percent_done = 100*(iteration/float(total))
    percent_str = ("{0:." + str(decimals) + "f}").format(percent_done)
    
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    
    print(f'\r{prefix} |{bar}| {percent_str}% {suffix}', end = print_end)
    
    # Print New Line on Complete
    if iteration == total: 
        print()