import numpy as np

def uniqueness(list1, list2):
    """
    Find the elements that are belonging to a group in both list1 and list2,
    where a group is defines as 2 or more elements that share the same
    value. 
    """
    assert len(list1) == len(list2)
    if len(list1)==1:
        return [1]
    l1_u = np.unique(list1)
    l2_u = np.unique(list2)
    unique = np.zeros(len(list1))
    kk = 0
    for u_l1 in l1_u:
        list1 == u_l1
        UK1 = np.where(list1 == u_l1)[0]
        UL1 = [pp in UK1 for pp in range(len(list1))]
        for u_l2 in l2_u:
            UK2 = np.where(list2 == u_l2)[0]
            UL2 = [pp in UK2 for pp in range(len(list1))]
            UUL = [UL1[pp]*UL2[pp] for pp in range(len(list1))]
            if len(np.where(np.array(UUL) != 0)[0]) == 0:
                continue 
            kk += 1
            unique += [kk*UUL[pp] for pp in range(len(list1))]
    # fill out zeros
    umax = np.max(unique)
    zeros = np.where(unique==0)[0]
    for ppk, pp in enumerate(zeros):
        unique[pp] += ppk+umax
    return unique.astype(int)

#add 'd0' to floating point number to avoid random trailing digits in Fortran input routines
def num2str(x):
    s = str(x)
    if s.find('e')<0:
        s += 'd0'
    return s
