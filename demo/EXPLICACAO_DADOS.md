# Explicação do Conteúdo dos Arquivos CSV

Cada um dos 5 arquivos CSV gerados para uma simulação (por exemplo, `Adsorption.csv`, `Heavy_Reflux.csv`, etc.) representa uma única etapa no processo PSA (Pressure Swing Adsorption). As linhas em cada arquivo correspondem a diferentes pontos no tempo durante essa etapa, e as colunas representam o estado do sistema em diferentes pontos espaciais dentro da coluna de adsorção.

## Estrutura das Colunas

A simulação utiliza um método de volumes finitos, o que significa que a coluna de adsorção é dividida em vários segmentos ou "nós". No script de demonstração, a coluna é dividida em **12 nós** (`N=10` mais dois nós de contorno). O nó 1 está no início da coluna (a entrada), e o nó 12 está no final (a saída).

Aqui está o que cada coluna nos arquivos CSV significa:

*   **Time**: Esta é a primeira coluna. Representa o tempo decorrido em segundos desde o início daquela etapa específica do processo.

*   **P_Node{i}**: (ex: `P_Node1`, `P_Node2`, ..., `P_Node12`)
    Esta é a **pressão adimensional** no nó `i`. Para obter a pressão real em Pascals (Pa), você precisa multiplicar este valor pela pressão de alimentação (`P_0`), que é uma das variáveis de processo definidas na simulação.

*   **y_Node{i}**: (ex: `y_Node1`, `y_Node2`, ..., `y_Node12`)
    Esta é a **fração molar de CO₂** na fase gasosa no nó `i`. Este valor é adimensional e varia de 0 a 1.

*   **x1_Node{i}**: (ex: `x1_Node1`, `x1_Node2`, ..., `x1_Node12`)
    Esta é a **quantidade adsorvida adimensional de CO₂** (componente 1) no material adsorvente no nó `i`. Representa a carga de CO₂ no adsorvente sólido.

*   **x2_Node{i}**: (ex: `x2_Node1`, `x2_Node2`, ..., `x2_Node12`)
    Esta é a **quantidade adsorvida adimensional de N₂** (componente 2) no material adsorvente no nó `i`.

*   **T_Node{i}**: (ex: `T_Node1`, `T_Node2`, ..., `T_Node12`)
    Esta é a **temperatura adimensional** no nó `i`. Para obter a temperatura real em Kelvin (K), você precisa multiplicar este valor pela temperatura de alimentação (`T_0`).

## Variáveis Adimensionais

Como mencionado, a maioria das variáveis é adimensional. O arquivo `README.md` do projeto fornece as definições exatas para convertê-las para suas unidades físicas:

*   **Pressão**: `P_real = P * P_0`
*   **Temperatura**: `T_real = T * T_0`
*   **Quantidade Adsorvida**: `q_i = x_i * q_s0` (onde `q_s0` é uma quantidade de referência adsorvida)

Analisando esses arquivos CSV, pode-se ver como a pressão, a temperatura e as concentrações de CO₂ e N₂ mudam ao longo do tempo e ao longo do comprimento da coluna de adsorção durante cada etapa do ciclo PSA.

