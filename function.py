def is_canonic_intron(intron):
    if len(intron) >= 4:
        return intron[:2].lower() == 'gt' and intron[-2:].lower() == 'ag'
    return False