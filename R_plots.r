## Dispersion plot -----------------------------------------------

# Dispersion plot is a scatterplot comparing mean expression level with 
# dispersions from mean to access the distribution. 
# Additional points show the fitted estimates, the final shrunk estimates and 
# outliers not used for fitting. It gives an idea of how well DESeq2 functioned.
# Note, the axes are in log scale.

plotDisp = function(.dds) {

    cols = mcols(.dds)

    select = cols$baseMean > 0

    df = data.frame(
        px = cols$baseMean[select],
        py = I(cols$dispGeneEst[select]),
        out = ifelse(cols$dispOutlier[select], "Outlier", "Estimate"),
        disp = I(dispersions(.dds)[select]),
        fit = cols$dispFit[select]
    )

    ymin = 10^floor(log10(min(df$py[df$py > 0], na.rm = TRUE)) - 0.1)

    library(ggplot2)

    ggplot(df, aes(px)) +
        geom_point(
            aes(
                y = py, col = out, alpha = out,
                shape = out 
            ),
            fill = "black", stroke = 1, size = 0.5
        ) +
        geom_point(
            aes(
                x = px, y = disp, col = "Final",
                alpha = "Estimate", shape = "Estimate"
            ),
            data = df[df$out == "Estimate", ], size = 0.5
        ) +
        geom_smooth(aes(y = fit, col = "Fitted"), size = 0.8, method = "gam") +
        scale_x_log10("Mean of normalized counts") +
        scale_y_log10("Dispersion") +
        scale_color_manual(
            "Values",
            values = c(
                "Estimate" = "black",
                "Outlier" = "#30a9f4",
                "Fitted" = "red",
                "Final" = "#17b978"
            )
        ) +
        scale_alpha_manual(
            values = c("Estimate" = .4, "Outlier" = 1),
            guide = F
        ) +
        scale_shape_manual(
            values = c("Estimate" = 16, "Outlier" = 21),
            guide = F
        ) +
        theme_classic() +
        guides(
            color = guide_legend(override.aes = list(
                linetype = c(0, 0, 1, 0),
                shape = c(20, 20, NA, 21),
                size = c(1.5, 1.5, 1, 2)
            ))
        )
}

## MA plot -----------------------------------------------

# Volcano plots compare log-fold change to mean expression level.
# Significant genes, at a provided alpha is identified.
# Note, the axes are in log scale.

plotMA = function(.data, .alpha, .ylim = c(2, -2)) {
    if(length(.ylim) == 1) .ylim = c(.ylim, -1 *.ylim)
    .data = na.omit(as.data.frame(.data))
    .data$sig = ifelse(.data$padj <= .alpha, "Yes", "No")
    .data$cex = ifelse(
        .data$log2FoldChange > .ylim[1], "r",
        ifelse(.data$log2FoldChange < .ylim[2], "l", "g")
    )
    .data$log2change = pmin(.ylim[1], pmax(.ylim[2], .data$log2FoldChange))
    library(ggplot2)
    ggplot(.data) +
        geom_point(
            aes(baseMean, log2change, fill = sig, shape = cex, col = sig),
            alpha = 0.4
        ) +
        scale_y_continuous("Log2 fold Change") +
        scale_x_log10("Mean normalized count") +
        scale_shape_manual(
            values = c("r" = 24, "l" = 25, "g" = 21),
            guide = F
        ) +
        scale_fill_manual(
            values = c("Yes" = "#D32F2F", "No" = "grey20"),
            guide = F
        ) +
        scale_color_manual(
            "Significant",
            breaks = c("Yes", "No"),
            values = c("Yes" = "#D32F2F", "No" = "grey20")
        )
}
```

## Volcano plot --------------------------------------------------

# Volcano plots compare log of p-value to log-fold change.
# Significant genes, at a provided alpha is identified.
# Note, the axes are in log scale.

plotVolcano = function(.data, .alpha, .xlim = 2, .ylim = 60) {
    .data = na.omit(as.data.frame(.data))
    .data$sig = .data$padj < .alpha & abs(.data$log2FoldChange) > 1
    .data$cex = ifelse( abs(.data$log2FoldChange) < .xlim, "i", "o")
    .data$log2change = pmax(-1 * .xlim, pmin(.xlim, .data$log2FoldChange))
    .data$pLog = pmin(.ylim, -1 * log10(.data$padj))
    .data$cex = ifelse(-1*log10(.data$padj) > .ylim, "u", .data$cex)
    library(ggplot2)
    ggplot(.data) +
        geom_point(aes(log2change, pLog, shape = cex, col = sig), alpha = 0.4) +
        scale_y_continuous(
            "adjusted p-value",
            breaks = 0:.ylim,
            expand = c(0, 0.1),
            labels = function(x) ifelse(x%%5 == 0, scales::math_format()(-1*x), "")
        ) +
        scale_x_continuous(
            "fold change",
            breaks = function(x) seq(ceiling(x[1]), x[2], 1),
            expand = c(0, 0.1),
            labels = function(x) ifelse(x%%2 == 0, scales::math_format(2^.x)(x), "")
        ) +
        scale_shape_manual(values = c("u" = 17, "o" = 18, "i" = 16), guide = F)+
        scale_color_manual(NULL,
            breaks = c(T, F),
            values = c("TRUE" = "#D32F2F", "FALSE" = "grey20"),
            labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
        ) +
        theme_classic(base_size = 22, base_family = "OpenSans") +
        theme(
            legend.position = "bottom",
            plot.title = element_text(hjust = 0.5),
            plot.background = element_rect(fill = "transparent", color = NA),
            legend.background = element_rect(fill = "transparent", color = NA),
            panel.background = element_rect(fill = "transparent", color = NA)
        )
}

## Individual gene count plot -------------------------------------

# Scatterplot showing raw counts of each replication for a individual gene
# across conditions.
# Note, the y-axis is in log scale.

countPlot = function(.dds, .gene, .cond) {
 library(ggplot2)
    library(DESeq2)

    count = counts(.dds)

    if (! .gene %in% rownames(count)) return(NULL)
    if (sum(count[.gene, ]) <= 20) return(NULL)

    d = do.call(rbind, lapply(.gene, function(x) {
        if (x %in% rownames(count)) {
            data.frame(count = counts(.dds)[x, ] + 1, gene = x)
        }
    }))

    df = cbind.data.frame(d, colData(.dds))

    p = ggplot(df, aes_string(.cond[1], "count")) +
        geom_point(position = position_jitter(w = 0.1, h = 0)) +
        scale_y_continuous(
            trans = "log2",
            breaks = scales::trans_breaks("log2", function(x) 2^x),
            labels = scales::trans_format("log2", scales::math_format(2^.x))
        ) +
        scale_x_discrete(NULL, labels = c(C = "Control", X = "Inoculated")) +
        theme_minimal(base_size = 15) +
        theme(
            panel.border = element_rect(size = 1, fill = "transparent"),
            strip.placement = "outside",
            plot.title = element_text(size = 16, hjust = 0.5, vjust = 1),
            plot.background = element_rect(color = "black", fill = NA),
            plot.margin = unit(c(2, 2,2 ,2), "mm"),
            legend.background = element_rect(fill = "transparent", color = NA),
            panel.background = element_rect(fill = "transparent", color = NA)
        )

    if (length(.cond) > 1) {
        p = p +
            facet_grid(reformulate(.cond[2]), labeller = as_labeller(
                c("2" = "2 dpi", "5" = "5 dpi")
            )) +
            stat_summary(
                aes_string(group = .cond[2]),
                fun = mean,
                geom = "line",
                col = "red"
            ) +
            scale_y_continuous(NULL,
                trans = "log2",
                # sec.axis = dup_axis(),
                breaks = scales::trans_breaks("log2", function(x) 2^x),
                labels = scales::trans_format("log2", scales::math_format(2^.x))
            )
    }

    p
}

## PCA plot ---------------------------------------------------

# 500 genes with highest variance will be used for generating PCA.
# Outliers samples are detected in the process.
# Dendrogram is produced for samples.

plotPCA = function(.dds, .cond= "condition", .ntop = 1000, .nclust = 2) {

    library(DESeq2)
    library(ggplot2)
    library(ggrepel)
    library(ggdendro)

    rld = rlog(.dds, blind = FALSE)
    assay = assay(.dds)

    # Select limited genes for PCA
    rv = rowVars(assay)
    select = order(rv, decreasing = T)[seq_len(min(.ntop, length(rv)))]

    # Compute PCA and variances
    pca = prcomp(t(assay[select, ]))
    pVar = pca$sdev^2 / sum(pca$sdev^2)

    U = pca$x
    print(Reduce(union,apply(U, 2, function(x) which( (abs(x - median(x)) / mad(x)) > 6 ))))
    
    eigen = data.frame(
        x = paste0("PC", 1:6),
        y = pVar[1:6] * 100
    )
    p1 = ggplot(eigen, aes(x, y)) + theme_void(base_size = 12) +
        geom_bar(stat = "identity", fill = "transparent", col = "black") +
        geom_text(aes(label = paste0(round(y, 2), "%")), nudge_y = 5) +
        ggtitle("Eigen values") +
        theme(
            plot.title = element_text(hjust = 0.5),
            plot.background = element_rect(color = "black", fill = "transparent")
        )

    # Identify clusters
    sdat = t(scale(t(assay[select, ])))
    tree = hclust(dist(t(sdat), method = "euclidean"))
    clust = cutree(tree, k = .nclust)
    tD = dendro_data(tree, type = "rectangle")
    p2 = ggplot() + theme_void(base_size = 12) +
        geom_segment(aes(x, y, xend = xend, yend = yend), segment(tD)) +
        geom_text(
            aes(x, y, label = label),
            label(tD), nudge_y = -18, angle = 90, hjust = 0
        ) +
        ggtitle("Dendrogram") +
        theme(
            plot.title = element_text(hjust = 0.5),
            plot.background = element_rect(color = "black", fill = "transparent")
        )

    out = apply(pca$x, 2, function(x) which(abs(x - median(x)) > (6 * mad(x))))
    if (length(out) > 1) print("Possible outlier found.")

    d = data.frame(
        PC1 = pca$x[, 1],
        PC2 = pca$x[, 2],
        name = colnames(.dds)
    )

    d$group = as.character(clust[d$name])

    colIndex = seq_len(length(unique(d$group)))
    col = c("#D32F2F", "#543dbd", "#a77820", "#19a8bb")[colIndex]

    sciLab = function(x) sapply(
        strsplit(scales::scientific_format()(x), "e\\+*0*"),
        function(i) ifelse(as.numeric(i[1]) == 0, "0", as.expression(
            bquote(.(i[1]) ~ "\u00D7" ~ 10^.(gsub("-0*", "-", i[2])))
        ))
    )

    p3 = ggplot(d, aes(PC1, PC2, col = group)) +
        geom_point(size = 3) +
        geom_text_repel(aes(label = as.character(name))) +
        scale_x_continuous(
            paste0("PC1: ", round(pVar[1] * 100), "% variance"),
            breaks = scales::pretty_breaks(n = 5),
            labels = sciLab
        ) +
        scale_y_continuous(
            paste0("PC2: ", round(pVar[2] * 100), "% variance"),
            breaks = scales::pretty_breaks(n = 5),
            labels = sciLab
        ) +
        scale_color_manual("Samples", values = col, guide = F) +
        scale_fill_manual(values = adjustcolor(col, 0.4), guide = F) +
        geom_hline(yintercept = 0, linetype = "dashed", col = "grey30") +
        geom_vline(xintercept = 0, linetype = "dashed", col = "grey30") +
        coord_fixed() +
        theme_classic(base_size = 18) +
        ggforce::geom_mark_ellipse(
            aes(fill = group, label = group, col = group),
            label.margin = margin(0, 0, 0, 0),
            label.buffer = unit(0, "mm"),
            con.type = "none"
        ) +
        theme(
            axis.line = element_blank(),
            text = element_text(family = "OpenSans"),
            panel.border = element_rect(fill = "transparent", size = 1),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA)
        )

    p0 = ggplot() + theme_void()
    cowplot::plot_grid(
        cowplot::plot_grid(
            p0, p1, p2, p0, nrow = 1, rel_widths = c(.45, .25, .25, .05)
        ), p3, ncol = 1, rel_heights = c(0.30, 0.65)
    )
}

## Heatmaps -----------------------------------------------------

# Heatmaps is generate for top (highest p-value) 50 genes by default, 
# with at least 2 fold change. Only selected genes is plotted
# when genelist is provided.

plotHeat = function(.dds, .genes = NULL, .cond, .ntop = 25, tree = F) {

    library(DESeq2)
    library(ggdendro)
    library(ggplot2)
    library(cowplot)

    # Normalizing the counts
    rld = rlog(.dds, blind = F)

    # remove genes with any samples with 0 reads
    counts = counts(.dds)
    zeros = sapply(seq_len(nrow(counts)), function(x) all(counts[x, ] != 0))

    # Select limited genes for PCA
    res = na.omit(results(.dds)[zeros, ])

    if (!is.null(.genes)) {
        res = res[match(.genes, rownames(res), nomatch = 0), ]
    } else {
        res = res[abs(res$log2FoldChange) >= 1, ]
        res = res[order(res$padj), ]
        if (.ntop) res = res[seq_len(min(.ntop, nrow(res))), ]
    }

    # Calculate distance
    res2 = assay(rld)
    rMat = t(scale(t(res2[match(rownames(res), rownames(res2)), ])))

    rDF = reshape2::melt(
        data.frame(
            rMat,
            gene = factor(rownames(rMat), levels = rev(rownames(rMat)))
        ),
        id = "gene",
        var = "sample"
    )

    if ("condition" %in% .cond) {
        rDF$condition = paste("Condition:", substr(rDF$sample, 1, 1))
    }

    if ("timePoint" %in% .cond) {
        rDF$timePoint = paste(
            "Timepoint:",
            ifelse(as.numeric(substr(rDF$sample, 2, 2)) > 3, 5, 2),
            "dpi"
        )
    }

    rDF$gene = paste0(
        annConv(rDF$gene, "Gene", "Gene_Desc"), " (",
        annConv(rDF$gene, "Gene", "Gene_ID"), ")"
    )

    ggplot(rDF, aes(sample, gene)) +
        geom_tile(aes(fill = pmin(pmax(value, -5), 5))) +
        scale_fill_gradientn(
            stringr::str_wrap("log2 count z-score", width = 12),
            colors = c("#4be787", "#4be787", "grey10", "#f63131", "#f63131"),
            values = scales::rescale(c(-5, -1.5, 0, 1.5, 5)),
            breaks = seq(-2, 2, .5),
            limits = c(-5, 5)
        ) +
        scale_x_discrete("Samples", position = "top", expand = c(0, 0)) +
        scale_y_discrete(
            "Genes",
            position = "right",
            labels = scales::wrap_format(30)
        ) +
        theme_classic() +
        facet_grid(
            reformulate(.cond, "."),
            scales = "free_x"
        ) +
        theme(
            axis.line = element_blank(),
            strip.placement = "outside",
            strip.background = element_blank(),
            panel.spacing = unit(0, "mm"),
            legend.key.height = unit(0.8, "in")
        )
}

## Venn diagram plots ----------------------------------------------

# plotVenn2 has coordinates for venn-diagrams with 2 sets,
# whereas plotVenn4 has coordinates for 4 sets, of which sets or two
# are mutaully exclusive.

plotVenn2 = function(.genelist, title = "", .labels = c("A", "B")) {

    plCir = data.frame(x = c(1, 2), y = 0, r = 1, group = .labels)

    plText = data.frame(
        x = c(0.5, 1.5, 2.5),
        y = 0,
        label = c(sum(.genelist$aOnly), sum(.genelist$ab), sum(.genelist$bOnly))
    )
    plText$percent = paste0(round(plText$label / 22252 * 100, 1), "%")

    plLabel = data.frame(
        x = c(0.3, 2.7),
        y = c(0.9, 0.9),
        hjust = c(1, 0),
        vjust = 1,
        label = .labels
    )

    library(ggplot2)

    ggplot() +
        ggforce::geom_circle(
            aes(x0 = x, y0 = y, r = r, fill = group),
            data = plCir, alpha = 0.5, col = "transparent"
        ) +
        geom_text(
            aes(x = x, y = y, label = label),
            plText,
            size = 10
        ) +
        geom_text(
            aes(x = x, y = y, label = percent),
            plText,
            size = 6,
            position = position_nudge(x = 0, y = -0.25)
        ) +
        geom_text(
            aes(x = x, y = y, label = label, vjust = vjust, hjust = hjust),
            plLabel,
            size = 7
        ) +
        ggforce::geom_circle(
            aes(x0 = x, y0 = y, r = r),
            data = plCir, alpha = 0, col = "black", size = 1.2
        ) +
        ggtitle(title) +
        scale_fill_manual(values = c("#ce3131", "#3694e0"), guide = F) +
        coord_equal() +
        theme_void() +
        theme(
            plot.title = element_text(
                size = 20,
                hjust = 0.5,
                face = "bold",
                vjust = 21
            ),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA)
        )

}

plotVenn4 = function(.genelist, title = "", .labels = c("A", "B", "C", "D")) {

    plCir = data.frame(
        x = c(0, 0, -1, 1),
        y = c(1, -1, 0, 0),
        r = 1,
        group = .labels
    )

    plText = data.frame(
        x = c(0, 0, -1.4, 1.4, -.5, .5, -.5, .5),
        y = c(1.4, -1.4, 0, 0, .5, .5, -.5, -.5),
        label = c(
            sum(.genelist$aUpOnly),
            sum(.genelist$aDownOnly),
            sum(.genelist$bUpOnly),
            sum(.genelist$bDownOnly),
            sum(.genelist$bothUp),
            sum(.genelist$aUpbDown),
            sum(.genelist$bDownbUp),
            sum(.genelist$bothDown)
            )
    )
    plText$percent = paste0(round(plText$label / 22252 * 100, 1), "%")

    plLabel = data.frame(
        x = c(-.9, 1, -1.4, 1.3),
        y = c(1.7, -1.5, 1.05, 1.05),
        hjust = c(1, 0, 1, 0),
        vjust = c(0, 1, .5, .5),
        label = .labels
    )

    library(ggplot2)

    ggplot() +
        ggforce::geom_circle(
            aes(x0 = x, y0 = y, r = r, fill = group),
            data = plCir, alpha = 0.5, col = "transparent"
        ) +
        geom_text(
            aes(x = x, y = y, label = label),
            plText,
            size = 10,
            position = position_nudge(x = 0, y = 0.05)
        ) +
        geom_text(
            aes(x = x, y = y, label = percent),
            plText,
            size = 6,
            position = position_nudge(x = 0, y = -0.15)
        ) +
        geom_text(
            aes(x = x, y = y, label = label, vjust = vjust, hjust = hjust),
            plLabel,
            size = 7
        ) +
        ggforce::geom_circle(
            aes(x0 = x, y0 = y, r = r),
            data = plCir, alpha = 0, size = 1.2, col = "black"
        ) +
        ggtitle(title) +
        scale_fill_manual(
            values = c("#ce3131", "#3694e0", "#51DF39", "#FDBF23"),
            guide = F
        ) +
        coord_fixed() +
        theme_void() +
        theme(
            plot.title = element_text(
                size = 20,
                hjust = 0.5,
                face = "bold"
            ),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA)
        )
}

## GO/KEGG terms Barplot --------------------------------------------

enrichBarPlot = function(.data, .ont, ...) {
    # Load ggplot
    library(ggplot2)

    # Alpha
    alpha = ifelse(is.null(list(...)$alpha), 0.05, list(...)$alpha)

    # Labs
    onts = list(
        "BP" = c("GO over-representation test", "Biological process"),
        "MF" = c("GO over-representation test", "Molecular function"),
        "CC" = c("GO over-representation test", "Cellular component"),
        "KEGG" = c("KEGG enrichment analysis", "Pathways")
    )

    ggplot(
            as.data.frame(.data),
            aes(Count, reorder(Description, rev(p.adjust)), fill = p.adjust)
        ) +
        geom_bar(stat = "identity") +
        scale_fill_gradient(
            "Adjusted p-value",
            low = "#D32F2F",
            high = "#303F9F",
            breaks = scales::pretty_breaks(n = 5),
            limits = c(0, alpha),
            guide = guide_colorbar(reverse = T)
        ) +
        labs(
            title = paste0(onts[[.ont]][1], "(", onts[[.ont]][2], ")"),
            x = "Number of genes",
            y = onts[[.ont]][2]
        ) +
        scale_y_discrete(labels = scales::wrap_format(30)) +
        theme_classic() +
        theme(
            plot.title = element_text(hjust = 0.5)
        )
}

## GO/KEGG terms Dotplot ----------------------------------------------

enrichDotPlot = function(.data, .ont, ...) {
    # Load ggplot
    library(ggplot2)

    # Alpha
    alpha = ifelse(is.null(list(...)$alpha), 0.05, list(...)$alpha)

    # Labs
    onts = list(
        "BP" = c("GO over-representation test", "Biological process"),
        "MF" = c("GO over-representation test", "Molecular function"),
        "CC" = c("GO over-representation test", "Cellular component"),
        "KEGG" = c("KEGG enrichment analysis", "Pathways")
    )

    ggplot(
            as.data.frame(.data),
            aes(
                sapply(GeneRatio, function(x) eval(parse(text = x))),
                reorder(Description, Count),
                color = p.adjust,
                size = Count
            )
        ) +
        geom_point() +
        scale_color_gradient(
            stringr::str_wrap("Adjusted p-value", width = 10),
            low = "#D32F2F",
            high = "#303F9F",
            breaks = function(x) seq(x[1], x[2], 0.01),
            limits = c(0, alpha),
            guide = guide_colorbar(reverse = T)
        ) +
        labs(
            title = paste0(onts[[.ont]][1], "(", onts[[.ont]][2], ")"),
            x = "Gene Ratio",
            y = onts[[.ont]][2]
        ) +
        scale_y_discrete(labels = scales::wrap_format(30)) +
        scale_size(
            stringr::str_wrap("Number of genes", width = 10),
            guide = guide_legend(reverse = T)
        ) +
        theme_classic() +
        theme(
            plot.title = element_text(hjust = 0.5)
        )
}