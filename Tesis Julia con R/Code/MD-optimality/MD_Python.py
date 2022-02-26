import pandas as pd
import numpy as np
import os
import itertools

def MD_Python(X, y, Xcand, nMod, p_mod, fac_mod, nFDes, max_int, g, Iter, nStart, top):
    # 1
    n = len(y)
    fac = len(X.columns) - 1

    Si = list()
    Xi = list()
    betai = list()
    gammai = list()
    efectos = list()

    models = np.zeros([nMod, fac])
    for i in range(0, nMod):
        aux = fac_mod.iloc[[i]]
        filtro = [a for a in aux.iloc[0] if a != 0]
        filtro[:] = [filt - 1 for filt in filtro]
        filtro = list(map(int, filtro))
        models[i, filtro] = 1

    Xfac = X.iloc[:, X.columns != "blk"]
    Xc   = Xcand.iloc[:, Xcand.columns != "blk"]

    # 2
    if max_int > 1:
        comb = list(itertools.combinations(range(1,fac + 1), 2))
        mat = np.zeros((len(models), len(comb)))
        for j in range(0, len(comb)):
            fac1 = comb[j][0]
            fac2 = comb[j][1]
            mat[np.where(models[:, fac1 - 1] + models[:, fac2 - 1] == 2), j] = 1
            Xfac['V' + str(j)] = Xfac.iloc[:, fac1 - 1]*Xfac.iloc[:, fac2 - 1]
            Xc['V' + str(j)] = Xc.iloc[:, fac1 - 1]*Xc.iloc[:, fac2 - 1]

        models = np.concatenate((models,mat),axis = 1)

    # 3
    if max_int > 2:
        comb = list(itertools.combinations(range(1,fac + 1), 3))
        mat = np.zeros((len(models), len(comb)))
        for j in range(0, len(comb)):
            fac1 = comb[j][0]
            fac2 = comb[j][1]
            fac3 = comb[j][2]
            mat[np.where(models[:, fac1 - 1] + models[:, fac2 - 1] + models[:, fac3 - 1] == 3), j] = 1
            Xfac['W' + str(j)] = Xfac.iloc[:, fac1 - 1]*Xfac.iloc[:, fac2 - 1]*Xfac.iloc[:, fac3 - 1]
            Xc['W' + str(j)] = Xc.iloc[:, fac1 - 1]*Xc.iloc[:, fac2 - 1]*Xc.iloc[:, fac3 - 1]

        models = np.concatenate((models,mat),axis = 1)

    Xfac.insert(0, 'X_1' ,X.iloc[:, 0])
    Xc.insert(0, 'Xcand', Xcand.iloc[:, 0], True)
    models = np.insert(models, 0, [1]*nMod, axis = 1)

    # 4
    for i in range(0, nMod):
        efectos.append(np.where(models[i, :] == 1))
        tam = sum([len(listElem) for listElem in efectos[i]])
        aux = Xfac.iloc[:, efectos[i][0]]
        aux.insert(0, '1',  [1]*n)
        Xi.append(aux)

        mat = np.zeros((tam + 1, tam + 1))
        if len(mat) > 1:
            coord = np.array([list(range(1, len(mat))), list(range(1, len(mat)))])
            for k in range(0, len(mat) - 1):
                mat[coord[0, k], coord[1, k]] = 1

        gammai.append((1 / g**2)*mat)
        betai.append(np.linalg.inv(gammai[i] + (np.transpose(Xi[i]) @ Xi[i])) @ np.transpose(Xi[i]) @ y)
        Si.append(np.transpose(y - Xi[i].to_numpy() @ betai[i]) @ (y - Xi[i].to_numpy() @ betai[i]) + np.transpose(betai[i]) @ gammai[i] @betai[i])

    # 5
    def MDr(extra):
        nex = len(extra)
        y_gorro_estrella = list()
        V_estrella = list()

        for i in range(0, nMod):
            Xiestrella = Xc.iloc[extra, efectos[i][0]]
            Xiestrella.insert(0, '1',  [1]*nex)
            Xiestrella = Xiestrella.to_numpy()

            y_gorro_estrella.append(Xiestrella @ betai[i])
            V_estrella.append(np.diag(np.ones(nex)) + Xiestrella @ np.linalg.inv(gammai[i] + (np.transpose(Xi[i]) @ Xi[i])) @ np.transpose(Xiestrella))

        MD = 0
        m = range(0, nMod)

        for i in m:
            for j in [elem for elem in m if m[elem] != i]:
                MD = MD + p_mod[i]*p_mod[j]*(-nex +
                                sum(np.diag(np.linalg.inv(V_estrella[j]) @ V_estrella[i])) +
                                (n - 1)*np.transpose(y_gorro_estrella[i] - y_gorro_estrella[j]) @
                                np.linalg.inv(V_estrella[j]) @ (y_gorro_estrella[i] - y_gorro_estrella[j]) / Si[i])
        MD = MD * 0.5

        return MD

    # 6
    df_MD = pd.DataFrame(columns = ['DesignPoints', 'MD'])
    df_MD.loc[len(df_MD)] = [0, 0]

    for j in range(0, nStart):
        extra = np.random.choice(range(0,len(Xcand)), size=nFDes, replace=True, p=None)
        iter = 1
        last_out = 0
        last_in = 1

        while (last_out != last_in) & (iter < Iter):
            dp = " ".join([str(int) for int in list(np.sort(extra))])

            aux = any(df_MD.iloc[:, 0] == dp)
            if aux != False:
                break

            df_MD.loc[len(df_MD)] = [dp] + [str(round(MDr(extra), 2))]

            op1 = list()
            for i in range(0, nFDes):
                op1.append(MDr(extra[np.arange(len(extra)) != i]))

            index = np.where(op1 == max(op1))[0][0]
            last_out = extra[index]
            extra = np.delete(extra, index)

            op2 = list()
            for i in range(0, len(Xcand)):
                op2.append(MDr(np.append(extra, i)))

            last_in = np.where(op2 == max(op2))[0][0]
            extra = np.append(extra, last_in)

            iter = iter + 1

    # 7
    df_MD['MD'] = df_MD['MD'].astype(float)
    df_MD = df_MD.sort_values(["MD"], ascending=False)
    df_MD = df_MD.iloc[0:top, :]

    names = (["F" + str(int) for int in list(range(1, fac + 1))])
    names.insert(0, "blk")
    X.set_axis(names, axis=1)

    Xcand.set_axis(names, axis = 1)

    return(df_MD)


# In[7]:
