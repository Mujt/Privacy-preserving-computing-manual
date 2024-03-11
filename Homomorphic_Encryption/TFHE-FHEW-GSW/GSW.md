## 1. 概述   
   GSW是第三代全同态加密的开山之作，是后续FHEW、TFHE的基础。其主要采用矩阵近似特征向量原理以及LWE来构造一种近似全同态（LFHE），再通过bootstrapping（Genntry09）操作实现全同态。本篇文章主要讲述GSW是如何构造全同态加密；会具体到一些基础的细节，比如LWE问题的原理、全同态加密的一般构造过程、bootstrapping原理等。本篇文章默认阅读者已经有一些密码学、线性代数基础，否则可以先移步密码学数学基础知识。
## 2. LWE(Learning With Error)
   首先，为了便于理解LWE是什么，先引入Learning without error的概念：
   比如，Alice有一组私钥(private key)，为一个向量，同时用私钥生成如下图右所示的方程组，即为公钥(public key):   
<div align = center><img src="https://github.com/Mujt/Privacy-preserving-computing-manual/assets/65994826/c61a0632-3912-43a7-a6fd-ce84b907e128" width="600px" height="360px" alt="Learning Without Error">
   <div align = left>显然，这里私钥就是公钥的解。我们需要私钥不能够被外界轻易获取，但是事实上外界通过公钥解方程很容易就能得到私钥，这就是Learning Without Error，我们不能通过这种方法来进行对数据的加密。
   而此时如果我们对公钥对应的方程中加入一些随机的噪声数据，事情似乎就变得有意思起来，如下图所示：
   <div align = center><img src="https://github.com/Mujt/Privacy-preserving-computing-manual/assets/65994826/2c93d6fe-6aa2-457f-93eb-d49cccfb1347" width="450px" height="360px" alt="Learning With Error">
  <img src="https://github.com/Mujt/Privacy-preserving-computing-manual/assets/65994826/127faf15-5b48-446c-9264-b68ace95ab82" width="450px" height="360px" alt="Learning With Error">
   <div align = left>此时，如果把这个已经添加过噪声的方程组作为公钥交由外部去操作，那么外部似乎是不那么容易获取我们的秘密信息--私钥了。事实上，这个充满错误的线性方程组是不可能求出正确解的，太多的方程限制了变量的值。这便是LWE问题的一般原理，外部无法通过一个充满错误的系统中获取正确的信息。   
      
   我们现在知道了这个原理，那么如何利用这个原理来进行一些数据的加密解密呢？仍以上述内容为例，给出Regev加密的流程。  
   需要注意的是，在密码学领域，一般运算都是在有限域下进行的（伽罗瓦域），这样做有一个好处是将数据约束在一定范围内。比如上述线性方程中，所有数据都要模一个素数89，得到下图所示公钥方程：
   <div align = center><img src="https://github.com/Mujt/Privacy-preserving-computing-manual/assets/65994826/9f35314f-1226-40a6-b815-4ffd30a4f79a" width="450px" height="360px" alt="Learning With Error">
   <div align = left>首先，Bob获取到了公钥中的三条数据：   

$77x+7y+28z+23w=8(mod89)$

$21x+19y+30z+48w=39(mod89)$ 

$4x+24y+33z+38w=20(mod89)$
   
   之后Bob将这三个方程进行相加，得到方程 $13x+50y+2z+20w=68(mod89)$ ，如果Bob想发送一个二进制的1，那么Bob则将该方程的右侧加上45（即89/2的向上取整，模数q是公开的），得到 $13x+50y+2z+20w=24(mod89)$，此时再将数据发送给Alice。此时Alice该怎么解密呢？事实上只要代入私钥即可！

   代入私钥之后，Alice得到 $13x+50y+2z+20w=69$ ，与Bob发来加密的数据进行对比，可知69与24差异较大，远远超出噪声可能的范围（一般应用时都会约束噪声的范围），因此可以得知Bob送来的是信息1。当Bob送来的数据为0时，也是直接代入私钥进行对比，此时差异在一个被约束的范围内，因此可以认为是0.

   以上便是一种典型的利用LWE进行数据加密解密的方法，即把一个bit的值映射到一个有限域的两头。

   总结来说，一个LWE问题实例会随机生成一个比较大的矩阵**A**（对应上述中的没有错误的线性方程组），和一个不公开的私密向量**s**（对应上述线性方程组的解）。给定一个**A**以及带有误差的乘积(**As**+**e**)，求出未知向量**s**的问题，叫做搜索LWE问题(DLWE)。还有一种LWE问题叫做决策LWE问题(DLWE)，是分辨看到的一组矩阵与向量到底是一个LWE实例(**A**,**As**+**e**)还是随机生成的(**A**,**v** $\in$ $Z$<sub>q</sub>)。一个合理构造的SLWE与DLWE在格密码学中都被定义为困难的问题。

   在我们懂了这种LWE的原理之后，我们继续学习GSW中如何用LWE问题来构造一个LFHE（Leveled Fully Homomorphic Encryption）近似全同态。
## 3. LFHE(Leveled Fully Homomorphic Encryption)
   GSW中使用的方法是矩阵近似特征向量。我们在矩阵特征向量公式中引入噪声，有以下公式：
   $$C \cdot \vec{s} = \mu \cdot \vec{s} + \vec{e}$$
   之后我们进行体系的构造：
   ### 密钥生成
   KeyGen：我们随机生成一个私密向量
   
   ### 加密算法ENC