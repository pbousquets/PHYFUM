/*
 * SubstitutionModel.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.beagle.evomodel.substmodel;

import dr.inference.model.Model;

/**
 * <b>model of sequence substitution (rate matrix)</b>.
 * provides a convenient interface for the computation of transition probabilities
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: SubstitutionModel.java,v 1.13 2005/05/24 20:25:58 rambaut Exp $
 */
public interface SubstitutionModel extends SubstitutionProcess, Model {
    // Combines SubstitutionProcess and Model
}