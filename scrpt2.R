#
#    script R para realizar experimentos fatoriais 2^k sobre o valor 
#    da evapotranspiração da cultura de referência
#
#    Copyright (C) 2015 Leonardo Cavalcante do Prado <leolleo.comp@gmail.com>    
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])

#instala a biblioteca r2html se não estiver instalada
if (!is.installed("R2HTML"))
	install.packages("R2HTML")

library(R2HTML)

#dados de entrada
dados <- read.table("input_dados.txt", sep=" ")
LIN_LIMIT = dim(dados)[1]

#rotulos dos dados
names(dados) <- c("Estacao", "Mes", "Ano", "Tmin_m", "Tmin_M", "Tmax_m", "Tmax_M", "URmin_m", "URmin_M", "URmax_m", "URmax_M", "PA_m", "PA_M", "Vv_m", "Vv_M", "RS_m", "RS_M")

#valores de ETo calculados a partir do programa em C++

#"legenda" da tabela
captions <- cbind(dados$Estacao, dados$Mes, dados$Ano)

for (i in 1:dim(captions)[1])
	captions[i,1] <- paste(c("A", captions[i,1]), collapse = "")

#alguns ajustes
dados$PA_m <- dados$PA_m / 10.0
dados$PA_M <- dados$PA_M / 10.0
dados$RS_m <- dados$RS_m / 1000.0
dados$RS_M <- dados$RS_M / 1000.0

#array de valores mínimos
input <- cbind(dados$Tmin_m, dados$Tmax_m, dados$URmin_m, dados$URmax_m, dados$PA_m, dados$Vv_m, dados$RS_m)
Lows  <- input

#array de valore máximos
input <- cbind(dados$Tmin_M, dados$Tmax_M, dados$URmin_M, dados$URmax_M, dados$PA_M, dados$Vv_M, dados$RS_M)
Highs <- input

#nome dos fatores
nomes <- as.vector(t(read.table("input_nomes.txt")))

#tabela de máximos e mínimos
mM_table_creator <- function(k)
{
	table <- NULL

	for (i in 0:(k-1))
		table <- cbind(table, rep(c(rep(-1, 2**i), rep(1, 2**i)), 2**(k - (i + 1))))

	for (i in 2:k) {
		comb    <- combn(1:k, i)
		lines   <- dim(comb)[1]
		columns <- dim(comb)[2]

		for (j in 1:columns) {
			tmp <- 1

			for (p in 1:lines)
				tmp <- tmp * table[,comb[p, j]]

			table <- cbind(table, tmp, deparse.level = 0)
		}
	}

	return (table)
}

# cabeçalho da tabela de máximos e mínimos
mM_table_header <- function(nomes, k)
{
	toSend <- nomes

	for (i in 2:k) {
		comb    <- combn(k, i)
		lines   <- dim(comb)[1]
		columns <- dim(comb)[2]

		for (j in 1:columns) {
			tmp <- paste(nomes[comb[,j]], collapse = "")
			toSend <- c(toSend, tmp)
		}
	}

	return (toSend)
}

# cabeçalho da tabela resumida
mM_resumed_table_header <- function(nomes, k)
{
	toSend <- nomes

	for (i in 2:k)
		toSend <- c(toSend, paste(c("C(", k, ", ", i, ")"), collapse = ""))

	return (toSend)
}

#cálculo dos efeitos
effects <- function(table, Y, k)
{
	q <- NULL

	for (i in 1:dim(table)[2])
		q <- c(q, sum(table[,i] * Y))

	q <- q / (2**(k-1))

	return (q)
}

#cálculo das variações
variations <- function(q, k)
{
	return ((2**(k-2)) * (q * q))
}

#SST
SST <- function(var)
{
	return (sum(var))
}

#variações explicadas
explained_fractions <- function(var, sst)
{
	return (100 * var / sst)
}

#variações explicadas resumido
resume_explained_fractions <- function(expfat, k)
{
	lim_inf <- k + 1
	toSend  <- expfat[1:k]

	for (i in 2:k) {
		toSend  <- c(toSend, sum(expfat[lim_inf:(lim_inf + choose(k, i) - 1)]))
		lim_inf <- lim_inf + choose(k, i)
	}

	return (toSend)
}

#retorna as coordenadas
coordenadas <- function(estacao)
{
	if (estacao == 307) return (c(-9.3833  * pi / 180.0,  370.5))
	if (estacao == 329) return (c(-8.5036  * pi / 180.0,  342.0))
	if (estacao == 330) return (c(-8.1325  * pi / 180.0,  374.0))
	if (estacao == 331) return (c(-8.3647  * pi / 180.0,  235.0))
	if (estacao == 337) return (c(-9.2861  * pi / 180.0,  100.0))
	if (estacao == 343) return (c(-7.0708  * pi / 180.0,  233.0))
	if (estacao == 345) return (c(-9.0331  * pi / 180.0,  402.0))
	if (estacao == 351) return (c(-8.6103  * pi / 180.0,  329.0))
	if (estacao == 354) return (c(-6.9742  * pi / 180.0,  156.0))
	if (estacao == 423) return (c(-9.6189  * pi / 180.0,  401.0))
	if (estacao == 424) return (c(-11.3289 * pi / 180.0,  755.0))
	if (estacao == 428) return (c(-10.4442 * pi / 180.0,  548.0))
	if (estacao == 435) return (c(-9.8336  * pi / 180.0,  453.0))
	if (estacao == 436) return (c(-10.9847 * pi / 180.0,  315.0))
	if (estacao == 440) return (c(-11.2050 * pi / 180.0,  453.0))
	if (estacao == 442) return (c(-10.5367 * pi / 180.0,  432.0))
	if (estacao == 443) return (c(-10.4553 * pi / 180.0,  637.0))
}

#juliano
jDay <- function(m, a, d = 15)
{
	j <- floor((275 * m / 9 - 30 + d) - 2)

	if (m < 3)
                j <- j + 2
	if (((a %% 400 == 0) || (a %% 4 == 0)) && ((a %% 100 != 0) && (m > 2)))
                j <- j + 1

	return (j)

}

#eo
eo <- function(T)
{
	return (0.6108 * exp((17.27 * T) / (T + 237.3)))
}

#ETo - Penman Monteith
ETo <- function (j, fi, z, Tmin, Tmax, URmin, URmax, PA, Vv, RS)
{
	u2 <- 0.0
	zm <- 10.0
	ds <- 0.409 * sin(2.0 * pi * j / 365.0 - 1.39)
	ws <- acos(-tan(fi) * tan(ds))
	dr <- 1.0 + 0.033 * cos(2.0 * pi * j / 365.0)
	Ra <- 37.58603 * dr * (ws * sin(fi) * sin(ds) + cos(fi) * cos(ds) * sin(ws))
	T  <- (Tmin + Tmax) / 2.0
	
	if (Vv > 0.08) # <?>
		u2 <- Vv * 4.87 / log(67.8 * zm - 5.42)

	es  = (eo(Tmin) + eo(Tmax)) / 2.0
	ea  = (eo(Tmin) * URmax + eo(Tmax) * URmin) / 200.0
	d   = 4098.0 * eo(T) / ((T + 237.3)**2)
	g   = 0.665E-3 * PA
	Rso = (0.75 + 2.0E-5 * z) * Ra
	Rnl = 4.903E-9 * (((Tmax + 273.16)**4) + ((Tmin + 273.16)**4)) / 2.0 * (0.34 - 0.14 * sqrt(ea)) * (1.35 * RS / Rso - 0.35)
	Rns = 0.77 * RS
	Rn  = Rns - Rnl

	return ((0.408 * d * Rn + g * 900.0 / (T + 273.0) * u2 * (es - ea)) / (d + g * (1.0 + 0.34 * u2)))
}

#dados do experimento
k               <- 7
maxMinTable     <- mM_table_creator(k)
header_1        <- c("Estacao", "Mes", "Ano", mM_table_header(nomes, k))
header_2        <- c("Estacao", "Mes", "Ano", mM_resumed_table_header(nomes, k))
output_stream_1 <- NULL
output_stream_2 <- NULL
resposta        <- NULL

#realizando o experimento para todas as estações
for (i in 1:dim(dados)[1]) {
	j     <- jDay(dados$Mes[i], dados$Ano[i])
	coord <- coordenadas(dados$Estacao[i])
        Y     <- NULL

	for (p in 1:dim(maxMinTable)[1]) {
                dadosEto <- NULL

                for (c in 1:k)
        		if (maxMinTable[p, c] > 0)
        			dadosEto  <- c(dadosEto, Highs[i, c])
                        else
                                dadosEto  <- c(dadosEto, Lows [i, c])

                Y <- c(Y, ETo(j, coord[1], coord[2], dadosEto[1], dadosEto[2], dadosEto[3], dadosEto[4], dadosEto[5], dadosEto[6], dadosEto[7]))
	}

	resposta <- rbind(resposta, Y)

	q       <- effects(maxMinTable, Y, k)
	vars    <- variations(q, k)
	expVar  <- explained_fractions(vars, SST(vars))
	expVarR <- resume_explained_fractions(expVar, k)

	output_stream_1 <- rbind(output_stream_1, expVar , deparse.level = 0)
	output_stream_2 <- rbind(output_stream_2, expVarR, deparse.level = 0)
}

write.table(resposta, file = "y_values", quote=FALSE, sep=";", row.names = FALSE)

#output_stream_1 = completo; output_stream_2 = resumido
output_stream_1 <- as.data.frame(output_stream_1, row.names = NULL)
output_stream_2 <- as.data.frame(output_stream_2, row.names = NULL)

#legenda
output_stream_1 <- cbind(captions, output_stream_1)
output_stream_2 <- cbind(captions, output_stream_2)

#Nomes das colunas
names(output_stream_1) <- header_1
names(output_stream_2) <- header_2

#escrevendo em .txt
write.table(output_stream_1, file = "output1.txt", quote = FALSE, sep = ";", row.names = FALSE)
write.table(output_stream_2, file = "output2.txt", quote = FALSE, sep = ";", row.names = FALSE)

#exportando pra html (somente o resumido)
HTMLInitFile(outdir = ".", filename = "output_htmlR", Title = "Resultados ETo", CSSFile = "")
HTML(output_stream_2, row.names = FALSE, nsmall = 3, classtable = "tablesorter", digits = 4, innerBorder = 1)
HTMLEndFile()
