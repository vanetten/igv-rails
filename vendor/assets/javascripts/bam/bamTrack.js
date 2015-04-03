/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2014 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/**
 * Created by turner on 2/24/14.
 */
var igv = (function (igv) {

    igv.BAMTrack = function (config) {

        igv.configTrack(this, config);

        this.visibilityWindow = config.visibilityWindow || 30000;     // 30kb default

        this.alignmentRowHeight = config.alignmentRowHeight || 14;

        this.coverageTrackHeight = config.coverageTrackHeight || 50;

        this.alignmentColor = config.alignmentColor || "rgb(185, 185, 185)";

        this.negStrandColor = config.negStrandColor || "rgba(150, 150, 230, 0.75)";
        this.posStrandColor = config.posStrandColor || "rgba(230, 150, 150, 0.75)";

        this.firstInfPairColor = "rgba(150, 150, 230, 0.75)";
        this.secondInPairColor = "rgba(230, 150, 150, 0.75)";

        this.deletionColor = config.deletionColor || "black";

        this.skippedColor = config.skippedColor || "rgb(150, 170, 170)";

        this.coverageColor = config.coverageColor || this.alignmentColor;

        this.maxHeight = config.maxHeight || 500;


        // divide the canvas into a coverage track region and an alignment track region
        this.alignmentRowYInset = 1;

        // alignment shading options
        //this.alignmentShading = config.alignmentShading || "none";
        this.alignmentShading = "none";

        // sort alignment rows
        this.sortOption = config.sortOption || { sort : "NUCLEOTIDE" };

        // filter alignments
        this.filterOption = config.filterOption || { name : "mappingQuality", params : [ 30, undefined ] };

        this.featureSource = new igv.BamSource(config);
    };

    igv.BAMTrack.alignmentShadingOptions = {

        none : function (bamTrack, alignment) {
            return bamTrack.alignmentColor;
        },

        strand : function (bamTrack, alignment) {
            return alignment.strand ? bamTrack.posStrandColor : bamTrack.negStrandColor;
        },

        firstOfPairStrand : function (bamTrack, alignment) {

            if (alignment.isPaired()) {

                if (alignment.isFistOfPair()) {
                    return alignment.strand ? bamTrack.posStrandColor : bamTrack.negStrandColor;
                }
                else if (alignment.isSecondOfPair()) {
                    return alignment.strand ? bamTrack.negStrandColor : bamTrack.posStrandColor;
                }
                else {
                    console.log("ERROR. Paired alignments are either first or second.")
                }

            } else {
                return bamTrack.alignmentColor;
            }

        }

    };

    igv.BAMTrack.counter = 1;

    igv.BAMTrack.filters = {

        noop : function () {
            return function (alignment) {
                return false;
            };
        },

        strand : function (strand) {
            return function (alignment) {
                return alignment.strand === strand;
            };
        },

        mappingQuality : function (lower, upper) {
            return function (alignment) {

                if (lower && alignment.mq < lower) {
                    return true;
                }

                if (upper && alignment.mq > upper) {
                    return true;
                }

                return false;
            }
        }
    };

    igv.BAMTrack.selectFilter = function (bamTrack, filterOption) {

        var a,
            b;

        if ("mappingQuality" === filterOption.name) {
            a = bamTrack.filterOption[ "params" ][ 0 ];
            b = bamTrack.filterOption[ "params" ][ 1 ];
            return igv.BAMTrack.filters[ filterOption.name ](a, b);
        }

        if ("strand" === filterOption.name) {
            a = bamTrack.filterOption[ "params" ][ 0 ];
            return igv.BAMTrack.filters[ filterOption.name ](a);
        }

        if ("noop" === filterOption.name) {
            return igv.BAMTrack.filters[ filterOption.name ]();
        }

        return undefined;
    };

    igv.BAMTrack.prototype.filterAlignments = function (filterOption, callback) {

        var pixelWidth,
            bpWidth,
            bpStart,
            bpEnd,
            filter;

        filter = igv.BAMTrack.selectFilter(this, filterOption);

        pixelWidth = 3 * this.trackView.canvas.width;
        bpWidth = Math.round(igv.browser.referenceFrame.toBP(pixelWidth));
        bpStart = Math.max(0, Math.round(igv.browser.referenceFrame.start - bpWidth / 3));
        bpEnd = bpStart + bpWidth;

        this.featureSource.getFeatures(igv.browser.referenceFrame.chr, bpStart, bpEnd, function (genomicInterval) {

            genomicInterval.packedAlignmentRows.forEach(function(alignmentRow){
                alignmentRow.alignments.forEach(function(alignment){
                    alignment.hidden = filter(alignment);
                });
            });

            callback();
        });
    };

    // Shift - Click to Filter alignments
    igv.BAMTrack.prototype.shiftClick = function (genomicLocation, event) {

        var myself = this;

        this.filterAlignments(this.filterOption, function () {
            myself.trackView.update();
            $(myself.trackView.viewportDiv).scrollTop(0);
        });

    };

    igv.BAMTrack.prototype.sortAlignmentRows = function (genomicLocation, sortOption) {

        var myself = this;

        this.featureSource.getFeatures(igv.browser.referenceFrame.chr, genomicLocation, (1 + genomicLocation), function (genomicInterval) {

            doSortAlignmentRows(genomicLocation, genomicInterval, sortOption);
            myself.trackView.update();
            $(myself.trackView.viewportDiv).scrollTop(0);
        });
    };

    function doSortAlignmentRows(genomicLocation, genomicInterval, sortOption) {

        var alignmentRows = genomicInterval.packedAlignmentRows,
            sequence = genomicInterval.sequence;

        if (sequence) {
            sequence = sequence.toUpperCase();
        } else {
            console.log("No sequence, no traversal. No discussion!");
            return;
        }

        alignmentRows.forEach(function(alignmentRow){
            alignmentRow.updateScore(genomicLocation, genomicInterval, sortOption);
        });

        alignmentRows.sort(function(a, b) {
            return a.score - b.score;
        });

    }

    // Alt - Click to Sort alignment rows
    igv.BAMTrack.prototype.altClick = function (genomicLocation, event) {

        this.sortAlignmentRows(genomicLocation, this.sortOption);

    };

    igv.BAMTrack.prototype.getFeatures = function (chr, bpStart, bpEnd, continuation, task) {

        // Don't try to draw alignments for windows > the visibility window
        if (igv.browser.trackViewportWidthBP() > this.visibilityWindow) {
            continuation({exceedsVisibilityWindow: true});
            return;
        }

        this.featureSource.getFeatures(chr, bpStart, bpEnd, continuation, task);

    };

    igv.BAMTrack.prototype.draw = function (options) {

        var genomicInterval = options.features,
            ctx = options.context,
            bpPerPixel = options.bpPerPixel,
            bpStart = options.bpStart,
            pixelWidth = options.pixelWidth,
            pixelHeight = options.pixelHeight,
            skippedColor = this.skippedColor,
            deletionColor = this.deletionColor,
            bpEnd = bpStart + pixelWidth * bpPerPixel + 1,
            myself = this,
            zoomInNoticeFontStyle = { font: '16px PT Sans', fillStyle: "rgba(64, 64, 64, 1)", strokeStyle: "rgba(64, 64, 64, 1)" };

        if (genomicInterval.exceedsVisibilityWindow) {

            for (var x = 200; x < pixelWidth; x += 400) {
                igv.Canvas.fillText.call(ctx, "Zoom in to see alignments", x, 20, zoomInNoticeFontStyle);
            }

            return;
        }

        if (genomicInterval) {
            drawCoverage(genomicInterval.coverageMap);
            drawAlignments(genomicInterval);
        }

        function drawCoverage(coverageMap) {
            var bp,
                x,
                y,
                w,
                h,
                refBase,
                i,
                len,
                item,
                accumulatedHeight,
                sequence;


            if (coverageMap.refSeq) sequence = coverageMap.refSeq.toUpperCase();

            // TODO -- why is covereageMap sometimes undefined !?
            if (coverageMap) {

                // paint backdrop color for all coverage buckets
                w = Math.max(1, 1.0 / bpPerPixel);
                for (i = 0, len = coverageMap.coverage.length; i < len; i++) {

                    bp = (coverageMap.bpStart + i);
                    if (bp < bpStart) continue;
                    if (bp > bpEnd) break;

                    item = coverageMap.coverage[i];
                    if (!item) continue;

                    h = (item.total / coverageMap.maximum) * myself.coverageTrackHeight;


                    y = myself.coverageTrackHeight - h;
                    x = (bp - bpStart) / bpPerPixel;

                    igv.Canvas.setProperties.call(ctx, {fillStyle: myself.coverageColor, trokeStyle: myself.coverageColor });
                    igv.Canvas.fillRect.call(ctx, x, y, w, h);

                    // coverage mismatch coloring
                    if (sequence) {

                        //if (171167156 === bp) {
                        //    console.log("bp " + igv.numberFormatter(bp));
                        //}

                        refBase = sequence[i];
                        if (item.isMismatch(refBase)) {

                            igv.Canvas.setProperties.call(ctx, { fillStyle: igv.nucleotideColors[ refBase ] });
                            igv.Canvas.fillRect.call(ctx, x, y, w, h);

                            accumulatedHeight = 0.0;
                            [ "A", "C", "T", "G" ].forEach(function(nucleotide){

                                var count,
                                    hh;

                                count = item[ "pos" + nucleotide] + item[ "neg" + nucleotide];


                                // non-logoritmic
                                hh = (count / coverageMap.maximum) * myself.coverageTrackHeight;

                                y = (myself.coverageTrackHeight - hh) - accumulatedHeight;
                                accumulatedHeight += hh;

                                igv.Canvas.setProperties.call(ctx, { fillStyle: igv.nucleotideColors[ nucleotide ] });
                                igv.Canvas.fillRect.call(ctx, x, y, w, hh);

                            });

                        }

                    }

                }

            }
        }

        function drawAlignments(genomicInterval) {

            var packedAlignmentRows = genomicInterval.packedAlignmentRows,
                sequence = genomicInterval.sequence;

            if (sequence) {
                sequence = sequence.toUpperCase();
            }


            // TODO -- how can packedAlignmentRows be undefined?
            if (packedAlignmentRows) {
                // alignment track
                packedAlignmentRows.forEach(function renderAlignmentRow(alignmentRow, i) {

                    var widthArrowHead = myself.alignmentRowHeight / 2.0,
                        yStrokedLine,
                        yRect,
                        height;

                    yRect = myself.alignmentRowYInset + myself.coverageTrackHeight + (myself.alignmentRowHeight * i) + 5;
                    height = myself.alignmentRowHeight - (2 * myself.alignmentRowYInset);
                    yStrokedLine = (height / 2.0) + yRect;

                    alignmentRow.alignments.forEach(function renderAlignment(alignment, indexAlignment) {

                        var xStart,
                            xEnd,
                            canvasColor;

                        if (true === alignment.hidden) {
                            return;
                        }

                        if ((alignment.start + alignment.lengthOnRef) < bpStart) return;
                        if (alignment.start > bpEnd) return;

                        xStart = (alignment.start - bpStart) / bpPerPixel;
                        xEnd = ((alignment.start + alignment.lengthOnRef) - bpStart) / bpPerPixel;

                        canvasColor = igv.BAMTrack.alignmentShadingOptions[ myself.alignmentShading ](myself, alignment);

                        if (alignment.blocks.length > 0) {

                            for (var c = 0; c < alignment.cigar.length; c++) {

                                if      ("D" === alignment.cigar.charAt( c )) {
                                    igv.Canvas.strokeLine.call(ctx, xStart, yStrokedLine, xEnd, yStrokedLine, {strokeStyle: deletionColor});
                                    break;
                                }
                                else if ("N" === alignment.cigar.charAt( c )) {
                                    igv.Canvas.strokeLine.call(ctx, xStart, yStrokedLine, xEnd, yStrokedLine, {strokeStyle: skippedColor});
                                    break;
                                }

                            }

                        }

                        igv.Canvas.setProperties.call(ctx, {   fillStyle: canvasColor , strokeStyle: canvasColor});

                        alignment.blocks.forEach(function (block, indexBlocks) {
                            var refOffset = block.start - bpStart,
                                seqOffset = block.start - genomicInterval.start,
                                xBlockStart = refOffset / bpPerPixel,
                                xBlockEnd = ((block.start + block.len) - bpStart) / bpPerPixel,
                                widthBlock = Math.max(1, xBlockEnd - xBlockStart),
                                blockSeq = block.seq.toUpperCase(),
                                refChar,
                                readChar,
                                readQual,
                                xBase,
                                widthBase,
                                colorBase,
                                x,
                                y;

                            if      ( true === alignment.strand && indexBlocks === alignment.blocks.length - 1) {

                                x = [
                                    xEnd,
                                    xEnd + widthArrowHead,
                                    xEnd,
                                    xEnd];

                                y = [
                                    yRect,
                                    yRect + (height/2.0),
                                    yRect + height,
                                    yRect];

                                igv.Canvas.fillPolygon.call(ctx, x, y, { fillStyle: canvasColor });
                            }
                            else if (false === alignment.strand && indexBlocks === 0) {

                                x = [
                                    xBlockStart,
                                    xBlockStart - widthArrowHead,
                                    xBlockStart,
                                    xBlockStart];

                                y = [
                                    yRect,
                                    yRect + (height/2.0),
                                    yRect + height,
                                    yRect];

                                igv.Canvas.fillPolygon.call(ctx, x, y, { fillStyle: canvasColor });
                            }

                            igv.Canvas.fillRect.call(ctx, xBlockStart, yRect, widthBlock, height, { fillStyle: "white" });
                            igv.Canvas.fillRect.call(ctx, xBlockStart, yRect, widthBlock, height, { fillStyle: canvasColor });

                            // Only do mismatch coloring if a refseq exists to do the comparison
                            if (sequence && blockSeq !== "*") {

                                for (var i = 0, len = blockSeq.length; i < len; i++) {

                                    readChar = blockSeq.charAt(i);
                                    refChar = sequence.charAt(seqOffset + i);
                                    if (readChar === "=") {
                                        readChar = refChar;
                                    }

                                    if (readChar === "X" || refChar !== readChar) {
                                        if (block.qual && block.qual.length > i) {
                                            readQual = block.qual[ i ];
                                            colorBase = shadedBaseColor(readQual, readChar, i + block.start);
                                        }
                                        else {
                                            colorBase = igv.nucleotideColors[readChar];
                                        }

                                        if (colorBase) {

                                            xBase = ((block.start + i) - bpStart) / bpPerPixel;
                                            widthBase = Math.max(1, 1 / bpPerPixel);
                                            igv.Canvas.fillRect.call(ctx, xBase, yRect, widthBase, height, { fillStyle: colorBase });
                                        }
                                    }
                                }
                            }

                        });
                    });
                });
            }
        }

    };

    igv.BAMTrack.prototype.popupData = function (genomicLocation, xOffset, yOffset) {

        var coverageMap = this.featureSource.genomicInterval.coverageMap,
            coverageMapIndex,
            coverage,
            packedAlignmentRows = this.featureSource.genomicInterval.packedAlignmentRows,
            packedAlignmentsIndex,
            alignmentRow,
            alignment,
            nameValues = [];

        packedAlignmentsIndex = Math.floor((yOffset - (this.alignmentRowYInset + this.coverageTrackHeight)) / this.alignmentRowHeight);

        if (packedAlignmentsIndex < 0) {

            coverageMapIndex = genomicLocation - coverageMap.bpStart;
            coverage = coverageMap.coverage[ coverageMapIndex ];

            if (coverage) {


                nameValues.push(igv.browser.referenceFrame.chr + ":" + igv.numberFormatter(1 + genomicLocation));

                nameValues.push({name: 'Total Count', value: coverage.total});

                // A
                tmp = coverage.posA + coverage.negA;
                if (tmp > 0)  tmp = tmp.toString() + " (" + Math.floor(((coverage.posA + coverage.negA) / coverage.total) * 100.0) + "%)";
                nameValues.push({name: 'A', value: tmp});


                // C
                tmp = coverage.posC + coverage.negC;
                if (tmp > 0)  tmp = tmp.toString() + " (" + Math.floor((tmp / coverage.total) * 100.0) + "%)";
                nameValues.push({name: 'C', value: tmp});

                // G
                tmp = coverage.posG + coverage.negG;
                if (tmp > 0)  tmp = tmp.toString() + " (" + Math.floor((tmp / coverage.total) * 100.0) + "%)";
                nameValues.push({name: 'G', value: tmp});

                // T
                tmp = coverage.posT + coverage.negT;
                if (tmp > 0)  tmp = tmp.toString() + " (" + Math.floor((tmp / coverage.total) * 100.0) + "%)";
                nameValues.push({name: 'T', value: tmp});

                // N
                tmp = coverage.posN + coverage.negN;
                if (tmp > 0)  tmp = tmp.toString() + " (" + Math.floor((tmp / coverage.total) * 100.0) + "%)";
                nameValues.push({name: 'N', value: tmp});

            }

        }

        else if (packedAlignmentsIndex < packedAlignmentRows.length) {

            alignmentRow = packedAlignmentRows[ packedAlignmentsIndex ];

            alignment = undefined;

            for (var i = 0, len = alignmentRow.alignments.length, tmp; i < len; i++) {

                tmp = alignmentRow.alignments[ i ];

                if (tmp.start <= genomicLocation && (tmp.start + tmp.lengthOnRef >= genomicLocation)) {
                    alignment = tmp;
                    break;
                }

            }

            if (alignment) {

                return alignment.popupData(genomicLocation);

            }
        }

        return nameValues;

    };

    igv.BAMTrack.prototype.popupMenuItems = function (popover) {

        var myself = this,
            menuItems = [],
            lut = { "none": "Color: None", "strand": "Color: Read Strand", "firstOfPairStrand": "Color: 1st of Pair Strand" },
            checkMark     = '<i class="fa fa-check fa-check-shim"></i>',
            checkMarkNone = '<i class="fa fa-check fa-check-shim fa-check-hidden"></i>',
            trackMenuItem = '<div class=\"igv-track-menu-item\">',
            trackMenuItemFirst = '<div class=\"igv-track-menu-item igv-track-menu-border-top\">';

        //menuItems.push(igv.colorPickerMenuItem(popover, this.trackView, "Set feature color", this.color));

        [ "none", "strand", "firstOfPairStrand" ].forEach(function(alignmentShading, index){

            var chosen,
                str;

            chosen = (0 === index) ? trackMenuItemFirst : trackMenuItem;
            str = (alignmentShading === myself.alignmentShading) ? chosen + checkMark + lut[ alignmentShading ] + '</div>' : chosen + checkMarkNone + lut[ alignmentShading ] + '</div>';

            menuItems.push({
                object: $(str),
                click: function () {
                    popover.hide();

                    myself.alignmentShading = alignmentShading;
                    myself.trackView.update();
                }
            });

        });

        return menuItems;

    };

    /**
     * Optional method to compute pixel height to accomodate the list of features.  The implementation below
     * has side effects (modifiying the samples hash).  This is unfortunate, but harmless.
     *
     * @param features
     * @returns {number}
     */
    igv.BAMTrack.prototype.computePixelHeight = function (features) {

        if (features.packedAlignmentRows) {
            return this.alignmentRowYInset + this.coverageTrackHeight + (this.alignmentRowHeight * features.packedAlignmentRows.length) + 5;
        }
        else {
            return this.height;
        }

    };

    function shadedBaseColor(qual, nucleotide, genomicLocation) {

        var color,
            alpha,
            minQ = 5,   //prefs.getAsInt(PreferenceManager.SAM_BASE_QUALITY_MIN),
            maxQ = 20,  //prefs.getAsInt(PreferenceManager.SAM_BASE_QUALITY_MAX);
            foregroundColor = igv.nucleotideColorComponents[nucleotide],
            backgroundColor = [255, 255, 255];   // White


        //if (171167156 === genomicLocation) {
        //    // NOTE: Add 1 when presenting genomic location
        //    console.log("shadedBaseColor - locus " + igv.numberFormatter(1 + genomicLocation) + " qual " + qual);
        //}

        if (!foregroundColor) return;

        if (qual < minQ) {
            alpha = 0.1;
        } else {
            alpha = Math.max(0.1, Math.min(1.0, 0.1 + 0.9 * (qual - minQ) / (maxQ - minQ)));
        }
        // Round alpha to nearest 0.1
        alpha = Math.round(alpha * 10) / 10.0;

        if (alpha >= 1) {
            color = igv.nucleotideColors[nucleotide];
        }
        else {
            color = "rgba(" + foregroundColor[0] + "," + foregroundColor[1] + "," + foregroundColor[2] + "," + alpha + ")";    //igv.getCompositeColor(backgroundColor, foregroundColor, alpha);
        }
        return color;
    }

    return igv;

})
(igv || {});
